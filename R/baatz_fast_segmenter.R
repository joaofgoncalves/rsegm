#' Fast Baatz-Schape multiresolution segmentation
#'
#' Performs a **fast** multiresolution image segmentation inspired by the
#' Baatz-Schape OBIA algorithm.
#'
#' Compared to \code{\link{baatz_segmenter}}, this function calls a streamlined
#' C++ implementation designed to be substantially faster on small to
#' medium rasters by using efficient region bookkeeping, a priority-queue of
#' merge candidates, and a disjoint-set union (DSU) structure to finalize labels.
#'
#' @param x A \code{\link[terra:SpatRaster-class]{SpatRaster}} or a character
#'   string giving the path to a raster readable by \pkg{terra}.
#'
#' @param scale_param Numeric scalar controlling the segmentation scale.
#'   Larger values tend to produce larger segments. Internally, the C++
#'   routine uses a threshold proportional to \code{scale_param^2} for merge
#'   acceptance.
#'
#' @param color_weight Numeric in \eqn{[0,1]} controlling the relative importance
#'   of spectral (color) heterogeneity versus shape. Higher values emphasize
#'   spectral similarity.
#'
#' @param compactness Numeric in \eqn{[0,1]} controlling the trade-off between
#'   compactness and smoothness in the shape criterion.
#'
#' @param band_weights Optional numeric vector with length equal to the number
#'   of bands in \code{x}. Provides per-band weights for the spectral component.
#'   If \code{NULL}, all bands are weighted equally.
#'
#' @param smooth Integer \eqn{\ge 0}. Optional mean filter window size (pixels)
#'   applied to the raster prior to mean-shift. If \code{0} (default), no
#'   pre-smoothing is applied. If \code{> 0}, a \code{smooth x smooth} mean
#'   filter is applied via \code{\link[terra:focal]{terra::focal()}}.
#'
#' @param output_file Optional character string. If provided, the resulting
#'   segmentation raster is written to this file via
#'   \code{\link[terra:writeRaster]{terra::writeRaster()}}.
#'
#' @param verbose Logical. If \code{TRUE}, prints progress messages from the
#'   R wrapper and enables verbose output from the underlying C++ routine.
#'   If \code{FALSE} (default), runs quietly.
#'
#' @details
#' The input raster is converted to a dense matrix of
#' \eqn{n_{pixels} \times n_{bands}}. Any missing values are replaced by the
#' corresponding band mean before segmentation.
#'
#' The C++ routine initializes one region per pixel and builds 4-neighbour
#' adjacencies, tracking the shared boundary length between neighbouring regions.
#' Candidate merges are scored with a *fusion factor* combining:
#' \itemize{
#'   \item a spectral term based on the increase in within-region variance
#'         (computed from per-region band sums and sums of squares), and
#'   \item a shape term that combines compactness and smoothness based on
#'         region perimeter and bounding-box perimeter, while accounting for
#'         the shared boundary length between the two regions.
#' }
#'
#' Merges are greedily applied in increasing fusion-factor order until no
#' candidate remains below the scale threshold. Final labels are remapped to
#' consecutive integers starting at 1.
#'
#' Label values are arbitrary and should be treated as categorical.
#'
#' @return A single-layer \code{SpatRaster} with integer segment labels in a
#'   layer named \code{"segments"}.
#'
#' @seealso \code{\link{baatz_segmenter}}, \code{\link{fh_segmenter}},
#'   \code{\link{fh_meanshift_segmenter}}
#'
#' @references
#' Baatz, M. & Schape, A. (2000).
#' \emph{Multiresolution Segmentation: An Optimization Approach for
#' High Quality Multi-Scale Image Segmentation}.
#' In: Strobl et al. (eds), Angewandte Geographische Informationsverarbeitung XII.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r <- terra::rast(sample_raster_path())
#'
#' seg <- baatz_fast_segmenter(
#'   r,
#'   scale_param = 50,
#'   color_weight = 0.9,
#'   compactness = 0.5,
#'   verbose = TRUE
#' )
#' plot(seg)
#' }
#'
#' @export
#'
baatz_fast_segmenter <- function(x,
                                 scale_param = 50,
                                 color_weight = 0.9,
                                 compactness = 0.5,
                                 band_weights = NULL,
                                 smooth=0L,
                                 output_file = NULL,
                                 verbose = FALSE) {

  vcat <- function(...) if (isTRUE(verbose)) cat(...)

  # ------------------------------------------------------------------
  # Input handling
  # ------------------------------------------------------------------
  if (is.character(x)) {
    vcat("Loading raster from: ", x, "\n")
    x <- terra::rast(x)
  }

  if (!inherits(x, "SpatRaster")) {
    stop("Input must be a SpatRaster object or file path", call. = FALSE)
  }

  nrows <- terra::nrow(x)
  ncols <- terra::ncol(x)
  nb    <- terra::nlyr(x)

  vcat("\nPreparing image data...\n")
  vcat(sprintf("  Dimensions: %d rows x %d cols x %d bands\n",
               nrows, ncols, nb))


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
  # Raster -> matrix (pixels x bands)
  # ------------------------------------------------------------------
  img_matrix <- terra::as.matrix(x, wide = FALSE)

  # ------------------------------------------------------------------
  # NA handling: replace with band means
  # ------------------------------------------------------------------
  if (anyNA(img_matrix)) {
    vcat("  Replacing NA values with band means\n")
    for (b in seq_len(nb)) {
      mu <- mean(img_matrix[, b], na.rm = TRUE)
      img_matrix[is.na(img_matrix[, b]), b] <- mu
    }
  }

  # ------------------------------------------------------------------
  # Run fast segmentation (C++ already supports verbosity)
  # ------------------------------------------------------------------
  vcat("\nRunning fast Baatz-Schape segmentation...\n")

  labels <- baatz_fast_segmenter_cpp(
    image              = img_matrix,
    nrows              = nrows,
    ncols              = ncols,
    scale_param        = scale_param,
    color_weight       = color_weight,
    compactness_weight = compactness,
    band_weights       = band_weights,
    verbose            = verbose
  )

  if (is.null(labels) || !is.atomic(labels)) {
    stop("baatz_fast_segmenter_cpp() returned no labels (NULL/non-vector).",
         call. = FALSE)
  }

  # ------------------------------------------------------------------
  # Labels -> raster
  # ------------------------------------------------------------------
  vcat("\nCreating output raster...\n")

  lab_mat <- matrix(labels, nrow = nrows, ncol = ncols, byrow = TRUE)
  seg <- terra::rast(
    lab_mat,
    crs    = terra::crs(x),
    extent = terra::ext(x)
  )
  names(seg) <- "segments"

  # ------------------------------------------------------------------
  # Optional write to disk
  # ------------------------------------------------------------------
  if (!is.null(output_file)) {
    vcat("Saving to: ", output_file, "\n")
    terra::writeRaster(seg, output_file, overwrite = TRUE)
  }

  if (isTRUE(verbose)) {
    nseg <- length(unique(labels))
    vcat("\nSegmentation complete!\n")
    vcat(sprintf("  Unique segments: %d\n", nseg))
    vcat(sprintf("  Average segment size: %.1f pixels\n",
                 (nrows * ncols) / max(1L, nseg)))
  }

  seg
}
