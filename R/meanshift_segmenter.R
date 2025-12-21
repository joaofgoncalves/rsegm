#' Mean-shift image segmentation
#'
#' Performs image segmentation using a **mean-shift filtering + clustering**
#' approach on multi-band rasters (RGB or multispectral).
#'
#' The method first applies mean-shift *range filtering* (edge-preserving
#' smoothing) in a local spatial neighborhood, then groups pixels into
#' connected components based on similarity of the filtered values. Small
#' components can be merged into the most similar neighboring component.
#'
#' This is a high-level R wrapper around a C++ implementation optimized for
#' raster grids and multi-band imagery.
#'
#' @param x A \code{\link[terra:SpatRaster-class]{SpatRaster}} with 1 or more
#'   bands.
#'
#' @param spatialr Integer \eqn{\ge 1}. Spatial radius (in pixels) of the local
#'   neighborhood used during mean-shift filtering. Larger values increase
#'   spatial smoothing and runtime.
#'
#' @param ranger Positive numeric. Range bandwidth controlling how strongly
#'   pixels are weighted by spectral similarity during mean-shift filtering.
#'   Smaller values preserve edges more aggressively; larger values smooth more.
#'
#' @param max_iter Integer \eqn{\ge 1}. Maximum number of mean-shift iterations
#'   per pixel.
#'
#' @param eps Positive numeric. Convergence tolerance for the mean-shift update.
#'   Iteration stops when the squared shift is \eqn{\le eps^2}.
#'
#' @param merge_thr Optional numeric threshold used during clustering of the
#'   filtered image. Neighboring pixels whose filtered spectral distance is
#'   \eqn{\le merge\_thr} are connected. If \code{NA} (default), the C++ routine
#'   uses \code{0.5 * ranger}.
#'
#' @param min_size Integer \eqn{\ge 1}. Minimum segment size (in pixels).
#'   Components smaller than \code{min_size} are merged into the most similar
#'   neighboring component (based on filtered spectral distance).
#'
#' @param eight Logical. If \code{TRUE}, uses 8-neighborhood connectivity for
#'   clustering and also includes diagonal offsets in the mean-shift spatial
#'   neighborhood. If \code{FALSE}, uses 4-neighborhood.
#'
#' @param scale_bands Logical. If \code{TRUE} (default), standardizes each band
#'   to z-scores \code{(x - mean) / sd} using \code{na.rm = TRUE}. This is often
#'   recommended for multispectral imagery with bands in different numeric
#'   ranges.
#'
#' @param smooth Integer \eqn{\ge 0}. Optional mean filter window size (pixels)
#'   applied to the raster prior to mean-shift. If \code{0} (default), no
#'   pre-smoothing is applied. If \code{> 0}, a \code{smooth x smooth} mean
#'   filter is applied via \code{\link[terra:focal]{terra::focal()}}.
#'
#' @param return_filtered Logical. If \code{FALSE} (default), returns only the
#'   segmentation labels raster. If \code{TRUE}, also returns the filtered
#'   multi-band raster produced by mean-shift.
#'
#' @param output_file Optional character string. If provided, the resulting
#'   segmentation raster is written to this file via
#'   \code{\link[terra:writeRaster]{terra::writeRaster()}}.
#'
#' @param verbose Logical. If \code{TRUE}, prints progress messages from the R
#'   wrapper. If \code{FALSE}, runs quietly and suppresses console output from
#'   the underlying C++ routine.
#'
#' @details
#' **Algorithm outline (C++ implementation):**
#'
#' \enumerate{
#' \item \strong{Validity mask:} A pixel is marked invalid if any band is
#' \code{NA} or not finite. Invalid pixels receive label \code{0} in C++ and
#' are converted to \code{NA} in the output raster.
#'
#' \item \strong{Mean-shift filtering:} For each valid pixel, the algorithm
#' iteratively updates a spectral vector \eqn{y} using a weighted mean of the
#' local neighborhood within \code{spatialr}. The weight is a product of:
#' \itemize{
#'   \item a Gaussian \emph{spatial} kernel (sigma \eqn{\approx spatialr/2}), and
#'   \item a Gaussian \emph{range} kernel based on spectral distance with
#'   bandwidth \code{ranger}.
#' }
#' Iteration stops after \code{max_iter} steps or when the update magnitude is
#' below \code{eps}.
#'
#' \item \strong{Clustering:} The filtered image is clustered by building
#' connected components on the grid. Neighboring pixels are linked when the
#' squared distance between their filtered spectral vectors is below
#' \code{merge_thr^2}. Connectivity follows \code{eight}.
#'
#' \item \strong{Minimum size enforcement:} Components smaller than
#' \code{min_size} are merged into the most similar neighboring component using
#' local neighbor search.
#'
#' \item \strong{Relabeling:} Components are relabeled to consecutive integers
#' \code{1..K}. Invalid pixels remain \code{NA} in the R output.
#' }
#'
#' Internally, pixel values are passed to C++ in a band-major layout
#' (all pixels of band 1, then all pixels of band 2, ...). The wrapper handles
#' the conversion from \code{terra} values.
#'
#' @return If \code{return_filtered = FALSE}, a single-layer
#' \code{\link[terra:SpatRaster-class]{SpatRaster}} of integer segment labels
#' named \code{"segment_id"}.
#'
#' If \code{return_filtered = TRUE}, a list with:
#' \describe{
#' \item{\code{segments}}{Single-layer labels \code{SpatRaster}.}
#' \item{\code{filtered}}{Multi-layer \code{SpatRaster} with the mean-shift
#'   filtered bands (same number of layers as \code{x}).}
#' }
#'
#' @seealso \code{\link{fh_segmenter}}, \code{\link{fh_meanshift_segmenter}},
#'   \code{\link{baatz_segmenter}}
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' r <- terra::rast(sample_raster_path())
#'
#' # MeanShift segmentation (labels only)
#' seg <- meanshift_segmenter(
#'   r,
#'   spatialr = 5,
#'   ranger = 1.0,
#'   min_size = 50,
#'   eight = TRUE,
#'   scale_bands = TRUE,
#'   smooth = 0,
#'   verbose = TRUE
#' )
#' plot(seg)
#'
#' # Also return the filtered image
#' res <- meanshift_segmenter(r, return_filtered = TRUE, verbose = FALSE)
#' plot(res$segments)
#' }
#'
#' @export
#'
meanshift_segmenter <- function(x,
                                spatialr = 5L,
                                ranger = 1.0,
                                max_iter = 10L,
                                eps = 1e-3,
                                merge_thr = NA_real_,
                                min_size = 30L,
                                eight = TRUE,
                                scale_bands = TRUE,
                                smooth = 0L,
                                return_filtered = FALSE,
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

  # ------------------------------------------------------------------
  # Optional smoothing
  # ------------------------------------------------------------------
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

  # ------------------------------------------------------------------
  # Raster -> values matrix (ncell x nb)
  # ------------------------------------------------------------------
  v <- terra::values(x, mat = TRUE)

  na_rows <- apply(v, 1L, function(row) any(is.na(row)))
  if (any(na_rows)) {
    vcat(sprintf("  Warning: %d pixels contain NA values\n", sum(na_rows)))
    vcat("  Note: NA pixels may yield NA segment ids in output\n")
  }

  # ------------------------------------------------------------------
  # Band scaling
  # ------------------------------------------------------------------
  if (isTRUE(scale_bands)) {
    vcat("  Scaling bands: z-score per band (NA-aware)\n")
    for (j in seq_len(ncol(v))) {
      mu  <- mean(v[, j], na.rm = TRUE)
      sdv <- stats::sd(v[, j], na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      v[, j] <- (v[, j] - mu) / sdv
    }
  }

  # IMPORTANT: band-major layout expected by C++
  img <- as.numeric(v)

  run_cpp <- function() {
    meanshift_segmenter_cpp(
      img_d = img,
      nrow  = nr,
      ncol  = nc,
      nb    = nb,
      spatialr = as.integer(spatialr),
      ranger   = ranger,
      max_iter = as.integer(max_iter),
      eps      = eps,
      merge_thr = merge_thr,
      min_size = as.integer(min_size),
      eight    = eight
    )
  }

  # ------------------------------------------------------------------
  # Run segmentation (silence C++ output if verbose = FALSE)
  # ------------------------------------------------------------------
  if (isTRUE(verbose)) {
    vcat("\nRunning MeanShift segmentation...\n")
    res <- run_cpp()
  } else {
    utils::capture.output(res <- run_cpp(), type = "output")
  }

  if (is.null(res) || is.null(res$labels)) {
    stop("meanshift_segmenter_cpp() returned no labels.", call. = FALSE)
  }

  # ------------------------------------------------------------------
  # Labels -> raster
  # ------------------------------------------------------------------
  vcat("\nCreating output raster...\n")

  seg <- terra::rast(x, nlyr = 1)
  lab <- res$labels
  lab[lab == 0] <- NA_integer_
  terra::values(seg) <- lab
  names(seg) <- "segment_id"

  if (!isTRUE(return_filtered)) {
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
      terra::writeRaster(seg, output_file, overwrite = TRUE)
    }
    return(seg)
  }

  # ------------------------------------------------------------------
  # Optional filtered raster stack
  # ------------------------------------------------------------------
  vcat("Creating filtered raster stack...\n")

  filt <- terra::rast(x)
  names(filt) <- names(x)
  terra::values(filt) <-
    t(matrix(res$filtered, nrow = nb, byrow = TRUE))  # ncell x nb

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
    terra::writeRaster(seg, output_file, overwrite = TRUE)
  }

  list(segments = seg, filtered = filt)
}

