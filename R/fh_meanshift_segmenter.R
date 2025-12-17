#' FH mean-shift hybrid geospatial image segmenter
#'
#' Performs high-quality, scalable image segmentation by combining a fast
#' graph-based region merging step (Felzenszwalb-Huttenlocher, FH) with a
#' region-level Mean-Shift refinement. This hybrid approach achieves
#' state-of-the-art segmentation quality while remaining computationally
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
                                   final_min_size = 80) {

  stopifnot(inherits(x, "SpatRaster"))

  if (smooth > 0) {
    w <- matrix(1, smooth, smooth)
    x <- terra::focal(x, w = w, fun = mean, na.policy = "omit", fillvalue = NA)
  }

  v <- terra::values(x, mat = TRUE)  # ncell x nb

  if (scale_bands) {
    for (j in seq_len(ncol(v))) {
      mu <- mean(v[, j], na.rm = TRUE)
      sdv <- stats::sd(v[, j], na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      v[, j] <- (v[, j] - mu) / sdv
    }
  }

  nr <- nrow(x)
  nc <- ncol(x)
  nb <- terra::nlyr(x)

  # IMPORTANT: band-major (all pixels of band1, then band2, ...)
  img <- as.numeric(v)

  lab <- fh_meanshift_segmenter_cpp(
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

  out <- terra::rast(x, nlyr = 1)
  terra::values(out) <- lab
  names(out) <- "segment_id"
  out
}

