#' FH mean-shift hybrid geospatial image segmenter
#'
#' Performs high-quality, scalable image segmentation by combining a fast
#' graph-based region merging step (Felzenszwalb-Huttenlocher, FH) with a
#' region-level Mean-Shift refinement. This hybrid approach achieves
#' good segmentation quality while remaining computationally
#' efficient for large, multi-band raster images.
#'
#' The algorithm proceeds in four main stages:
#' \enumerate{
#'   \item Optional spatial smoothing and per-band normalization.
#'   \item Initial over-segmentation using the FH graph-based method
#'         (fast, noise-robust, edge-preserving).
#'   \item Region-level Mean-Shift clustering in reduced spectral space
#'         to merge spectrally similar regions.
#'   \item Final region merging and minimum-size enforcement.
#' }
#'
#' Compared to pixel-level Mean-Shift, this hybrid strategy is substantially
#' faster and more memory-efficient, while producing spatially coherent,
#' object-like segments suitable for OBIA workflows.
#'
#' @param x SpatRaster.
#'   Input image to segment. May be multi-band (e.g., RGB, VNIR, VNIR+SWIR).
#'
#' @param scale_bands logical, default TRUE.
#'   If TRUE, each band is centered and scaled to unit variance prior to
#'   segmentation. Strongly recommended for multi-sensor or multi-band data.
#'
#' @param smooth integer, default 3.
#'   Size of a square mean-filter kernel applied before segmentation.
#'   Set to 0 to disable smoothing. Values between 3 and 5 often reduce
#'   speckle and improve boundary stability.
#'
#' @param fh_k numeric, default 0.5.
#'   Scale parameter for the FH region-merging step. Larger values produce
#'   fewer, larger initial regions; smaller values preserve fine detail.
#'
#' @param fh_min_size integer, default 20.
#'   Minimum allowed size (in pixels) for FH regions. Smaller regions are
#'   merged during the FH cleanup stage.
#'
#' @param eight logical, default TRUE.
#'   If TRUE, uses 8-neighborhood connectivity; otherwise 4-neighborhood.
#'
#' @param ms_dim integer, default 0.
#'   Number of dimensions/bands to use in mean shift step.
#'
#' @param ms_ranger numeric, default 0.15.
#'   Spectral (range) bandwidth for Mean-Shift clustering of regions.
#'   Smaller values preserve spectral distinctions; larger values encourage
#'   region merging.
#'
#' @param ms_hs numeric, default 12.
#'   Spatial bandwidth (in pixels) used during Mean-Shift refinement.
#'   Controls spatial coherence of merged regions.
#'
#' @param ms_max_iter integer, default 10.
#'   Maximum number of Mean-Shift iterations per region.
#'
#' @param ms_eps numeric, default 1e-3.
#'   Convergence tolerance for Mean-Shift mode estimation.
#'
#' @param mode_merge numeric, default 0.6.
#'   Threshold controlling the merging of nearby Mean-Shift modes.
#'   Higher values result in more aggressive merging of similar regions.
#'
#' @param final_min_size integer, default 80.
#'   Minimum segment size (in pixels) enforced after Mean-Shift refinement.
#'   Remaining smaller segments are merged into spectrally closest neighbors.
#'
#' @param output_file Optional character string. If provided, the resulting
#'   segmentation raster is written to this file via
#'   \code{\link[terra:writeRaster]{terra::writeRaster()}}.
#'
#' @param verbose Do progress messages? (default: TRUE)
#'
#' @return SpatRaster.
#'   A single-layer raster where each cell contains an integer segment ID.
#'
#' @details
#' This function is optimized for large images and is well suited for
#' object-based image analysis (OBIA) of satellite or UAV imagery.
#' For very large rasters (tens to hundreds of millions of pixels),
#' consider tiling.
#'
#' @seealso
#'   \code{\link{fh_segmenter}}
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' img <- rast("multispectral_image.tif")
#'
#' seg <- fh_meanshift_segmenter(
#'   img,
#'   fh_k = 0.4,
#'   fh_min_size = 30,
#'   ms_ranger = 0.12,
#'   ms_hs = 10,
#'   final_min_size = 100
#' )
#'
#' plot(seg)
#' }
#'
#' @export
fh_meanshift_segmenter <- function(x,
                                   scale_bands = TRUE,
                                   smooth = 3,
                                   fh_k = 0.5,
                                   fh_min_size = 20,
                                   eight = TRUE,
                                   ms_dim = 0,
                                   ms_ranger = 0.15,
                                   ms_hs = 12,
                                   ms_max_iter = 10,
                                   ms_eps = 1e-3,
                                   mode_merge = 0.6,
                                   final_min_size = 80,
                                   output_file=NULL,
                                   verbose = TRUE) {

  vcat <- function(...) if (isTRUE(verbose)) cat(...)

  if (!inherits(x, "SpatRaster")) {
    stop("Input must be a SpatRaster object", call. = FALSE)
  }

  nr <- terra::nrow(x)
  nc <- terra::ncol(x)
  nb <- terra::nlyr(x)

  vcat("\nPreparing image data...\n")
  vcat(sprintf("  Dimensions: %d rows x %d cols x %d bands\n", nr, nc, nb))

  # Optional smoothing (mean filter)
  if (smooth > 0) {
    vcat(sprintf("  Smoothing: mean filter (%dx%d)\n", smooth, smooth))
    w <- matrix(1, smooth, smooth)
    x <- terra::focal(
      x,
      w = w,
      fun = mean,
      na.policy = "omit",
      fillvalue = NA
    )
  }

  # Raster -> values matrix (ncell x nb)
  v <- terra::values(x, mat = TRUE)

  # NA info (the C++ can keep invalid as NA if keep_invalid_na=TRUE)
  na_rows <- apply(v, 1L, function(row) any(is.na(row)))
  if (any(na_rows)) {
    vcat(sprintf("  Warning: %d pixels contain NA values\n", sum(na_rows)))
    vcat("  Note: NA pixels may yield NA segment ids in output\n")
  }

  # Band scaling
  if (isTRUE(scale_bands)) {
    vcat("  Scaling bands: z-score per band (NA-aware)\n")
    for (j in seq_len(ncol(v))) {
      mu  <- mean(v[, j], na.rm = TRUE)
      sdv <- stats::sd(v[, j], na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      v[, j] <- (v[, j] - mu) / sdv
    }
  }

  # IMPORTANT: band-major (all pixels of band1, then band2, ...)
  img <- as.numeric(v)

  run_cpp <- function() {
    fh_meanshift_segmenter_cpp(
      img_d = img,
      nrow  = nr,
      ncol  = nc,
      nb    = nb,
      fh_k = fh_k,
      fh_min_size = as.integer(fh_min_size),
      eight = eight,
      ms_dim = as.integer(ms_dim),
      ms_ranger = ms_ranger,
      ms_hs = ms_hs,
      ms_max_iter = as.integer(ms_max_iter),
      ms_eps = ms_eps,
      mode_merge = mode_merge,
      final_min_size = as.integer(final_min_size),
      keep_invalid_na = TRUE
    )
  }

  # Run segmentation (silence C++ output when verbose = FALSE)
  if (isTRUE(verbose)) {
    vcat("\nRunning FH + MeanShift segmentation...\n")
    lab <- run_cpp()
  } else {
    utils::capture.output(lab <- run_cpp(), type = "output")
  }

  if (is.null(lab) || !is.atomic(lab)) {
    stop("fh_meanshift_segmenter_cpp() returned no labels (NULL/non-vector).", call. = FALSE)
  }

  vcat("\nCreating output raster...\n")
  out <- terra::rast(x, nlyr = 1)
  terra::values(out) <- lab
  names(out) <- "segment_id"

  if (isTRUE(verbose)) {
    nseg <- length(unique(stats::na.omit(lab)))
    vcat("\nSegmentation complete!\n")
    vcat(sprintf("  Unique segments: %d\n", nseg))
    vcat(sprintf("  Average segment size: %.1f pixels\n",
                 (nr * nc) / max(1L, nseg)))
  }

  # ------------------------------------------------------------------
  # optional write to disk
  # ------------------------------------------------------------------
  if (!is.null(output_file)) {
    vcat("Saving to: ", output_file, "\n")
    terra::writeRaster(lab, output_file, overwrite = TRUE)
  }

  out
}
