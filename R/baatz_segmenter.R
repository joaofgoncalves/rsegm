
#' Baatz-Schape multiresolution segmentation
#'
#' Performs multiresolution image segmentation following the
#' Baatz-Schape algorithm, commonly used in object-based image analysis (OBIA).
#' The algorithm iteratively merges adjacent pixels or regions based on
#' spectral heterogeneity and shape constraints.
#'
#' This function is a high-level R wrapper around a lower-level
#' Baatz-Schape segmentation implementation, handling raster I/O,
#' preprocessing, and conversion of the output back to a `SpatRaster`.
#'
#' @param x A `SpatRaster` object or a character string giving the path
#'   to a raster file readable by \pkg{terra}.
#'
#' @param scale_param Numeric scalar controlling the segmentation scale.
#'   Larger values result in larger segments. Roughly corresponds to the
#'   maximum allowed increase in heterogeneity when merging regions.
#'
#' @param color_weight Numeric value in \eqn{[0,1]} giving the relative
#'   importance of spectral (color) heterogeneity versus shape.
#'   Higher values emphasize spectral similarity.
#'
#' @param compactness Numeric value in \eqn{[0,1]} controlling the trade-off
#'   between compactness and smoothness in the shape criterion.
#'
#' @param band_weights Optional numeric vector of length equal to the
#'   number of bands in `x`, giving per-band weights for the spectral
#'   heterogeneity term. If `NULL`, all bands are weighted equally.
#'
#' @param smooth Integer \eqn{\ge 0}. Optional mean filter window size (pixels)
#'   applied to the raster prior to mean-shift. If \code{0} (default), no
#'   pre-smoothing is applied. If \code{> 0}, a \code{smooth x smooth} mean
#'   filter is applied via \code{\link[terra:focal]{terra::focal()}}.
#'
#' @param output_file Optional character string. If provided, the resulting
#'   segmentation raster is written to this file using
#'   \code{\link[terra]{writeRaster}}.
#'
#' @param verbose Do progress messages? (default: TRUE)
#'
#' @details
#' Input rasters are internally converted to a matrix of
#' \eqn{n_{pixels} \times n_{bands}}. Pixels containing `NA` values in any
#' band are replaced by the corresponding band mean prior to segmentation.
#'
#' The segmentation labels are returned as integer region identifiers.
#' Label values are arbitrary and should be treated as categorical.
#'
#' This implementation is intended for small to medium-sized rasters.
#' For very large images, consider tiling strategies or out-of-core
#' processing.
#'
#' @return A single-layer `SpatRaster` with integer segment labels.
#'
#' @references
#' Baatz, M. & Schape, A. (2000).
#' \emph{Multiresolution Segmentation: An Optimization Approach for
#' High Quality Multi-Scale Image Segmentation}.
#' In: Strobl et al. (eds), Angewandte Geographische Informationsverarbeitung XII.
#'
#' @seealso
#' \code{\link{fh_segmenter}},
#' \code{\link{fh_meanshift_segmenter}},
#' \code{\link[terra]{rast}}
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r <- rast(system.file("extdata", "sample_raster.tif", package = "rsegm"))
#'
#' seg <- baatz_segmenter(
#'   r,
#'   scale_param = 50,
#'   color_weight = 0.9,
#'   compactness = 0.5
#' )
#'
#' plot(seg)
#' }
#'
#' @export
#'
baatz_segmenter <- function(x,
                            scale_param = 50,
                            color_weight = 0.9,
                            compactness = 0.5,
                            band_weights = NULL,
                            smooth=0L,
                            output_file = NULL,
                            verbose = TRUE) {

  # ------------------------------------------------------------------
  # helper for conditional messages
  # ------------------------------------------------------------------
  vcat <- function(...) {
    if (isTRUE(verbose)) {
      cat(...)
    }
  }

  # ------------------------------------------------------------------
  # input handling
  # ------------------------------------------------------------------
  if (is.character(x)) {
    vcat("Loading raster from: ", x, "\n")
    x <- terra::rast(x)
  }

  if (!inherits(x, "SpatRaster")) {
    stop("Input must be a SpatRaster object or file path", call. = FALSE)
  }

  # dimensions
  nrows  <- terra::nrow(x)
  ncols  <- terra::ncol(x)
  nbands <- terra::nlyr(x)

  vcat("\nPreparing image data...\n")
  vcat(sprintf("  Dimensions: %d rows x %d cols x %d bands\n",
               nrows, ncols, nbands))


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
  # raster matrix
  # ------------------------------------------------------------------
  img_matrix <- terra::as.matrix(x, wide = FALSE)

  # ------------------------------------------------------------------
  # NA handling
  # ------------------------------------------------------------------
  na_mask <- apply(img_matrix, 1L, function(row) any(is.na(row)))

  if (any(na_mask)) {
    vcat(sprintf("  Warning: %d pixels contain NA values\n", sum(na_mask)))
    vcat("  Replacing NA with band means\n")

    for (b in seq_len(nbands)) {
      mu <- mean(img_matrix[, b], na.rm = TRUE)
      img_matrix[is.na(img_matrix[, b]), b] <- mu
    }
  }

  # ------------------------------------------------------------------
  # segmentation (silence C++ output if needed)
  # ------------------------------------------------------------------
  run_cpp <- function() {
    baatz_segmenter_cpp(
      image               = img_matrix,
      nrows               = nrows,
      ncols               = ncols,
      scale_param         = scale_param,
      color_weight        = color_weight,
      compactness_weight  = compactness,
      band_weights        = band_weights
    )
  }

  if (isTRUE(verbose)) {
    vcat("\nRunning Baatz-Schape segmentation...\n")
    labels <- run_cpp()
  } else {
    utils::capture.output(labels <- run_cpp(), type = "output")
    #labels <- attr(labels, "value")
  }

  # ------------------------------------------------------------------
  # labels the raster
  # ------------------------------------------------------------------
  vcat("\nCreating output raster...\n")

  label_matrix <- matrix(labels,
                         nrow  = nrows,
                         ncol  = ncols,
                         byrow = FALSE)

  seg_rast <- terra::rast(
    t(label_matrix),
    crs    = terra::crs(x),
    extent = terra::ext(x)
  )

  names(seg_rast) <- "segments"

  # ------------------------------------------------------------------
  # optional write to disk
  # ------------------------------------------------------------------
  if (!is.null(output_file)) {
    vcat("Saving to: ", output_file, "\n")
    terra::writeRaster(seg_rast, output_file, overwrite = TRUE)
  }

  # ------------------------------------------------------------------
  # summary
  # ------------------------------------------------------------------
  if (isTRUE(verbose)) {
    nseg <- length(unique(labels))
    vcat("\nSegmentation complete!\n")
    vcat(sprintf("  Unique segments: %d\n", nseg))
    vcat(sprintf("  Average segment size: %.1f pixels\n",
                 (nrows * ncols) / nseg))
  }

  seg_rast
}


