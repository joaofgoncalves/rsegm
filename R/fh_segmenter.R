
#' Felzenszwalb-Huttenlocher (FH) graph-based image segmentation
#'
#' Segments a multi-band image using the **Felzenszwalb-Huttenlocher** graph-based
#' region merging algorithm (often used as an efficient OBIA-style "superpixel"
#' generator). Each pixel is a node in a grid graph; edges connect neighboring
#' pixels (4- or 8-neighborhood) with weights given by spectral distance.
#' Regions are merged in increasing edge-weight order subject to an adaptive
#' internal-difference criterion controlled by \code{k}, with an optional cleanup
#' pass that enforces a minimum region size.
#'
#' This wrapper accepts a \code{\link[terra]{SpatRaster}} and returns a single-layer
#' label raster with integer segment IDs. It optionally applies a mean (box) filter
#' prior to segmentation and z-score scales each band to stabilize the distance
#' threshold across bands and sensors.
#'
#' @details
#' **Algorithm behavior**
#' \itemize{
#'   \item \strong{\code{k}} controls the tendency to merge: larger values typically produce
#'   larger segments.
#'   \item \strong{\code{min_size}} enforces a minimum mapping unit by merging small regions
#'   to the most similar neighbor in a post-processing step.
#'   \item \strong{\code{eight}} toggles 8-neighborhood (Queen) vs 4-neighborhood (Rook)
#'   connectivity; 8-neighborhood usually yields more "diagonal-safe" objects.
#'   \item \strong{\code{scale_bands}} applies per-band z-score scaling using \code{na.rm=TRUE}.
#' }
#'
#' **NA handling**
#' Pixels with \code{NA} in any band are treated as invalid by the underlying C++
#' implementation and returned as \code{NA} in the output label raster.
#' If your scene contains large \code{NA} areas, consider masking/cropping/filling
#' prior to segmentation to avoid creating many invalid edge cases.
#'
#' **Performance notes**
#' The heavy lifting is implemented in Rcpp/C++ and is designed for large rasters.
#' For very large scenes, prefer tiling.
#'
#' @param x A \code{\link[terra]{SpatRaster}} with one or more layers (bands).
#'   Bands should be numeric and represent co-registered imagery (e.g., RGB, VNIR/SWIR).
#' @param k Numeric scalar. Scale/threshold parameter controlling region merging.
#'   Larger \code{k} generally yields fewer, larger segments.
#' @param min_size Integer scalar. Minimum allowed segment size (in pixels).
#'   Segments smaller than this are merged during a cleanup pass.
#' @param eight Logical. If \code{TRUE}, use 8-neighborhood connectivity; if \code{FALSE},
#'   use 4-neighborhood connectivity.
#' @param scale_bands Logical. If \code{TRUE} (recommended), each band is z-score scaled
#'   using its mean and standard deviation computed with \code{na.rm=TRUE}.
#' @param smooth Integer scalar. If \code{> 0}, applies a mean (box) filter of size
#'   \code{smooth x smooth} to each band before segmentation. Use small odd values
#'   (e.g., 3 or 5). Set to \code{0} to disable smoothing.
#'
#' @return A single-layer \code{\link[terra]{SpatRaster}} with integer segment IDs
#'   in the layer \code{"segment_id"}. Invalid pixels are \code{NA}.
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[terra]{rast}}, \code{\link[terra]{values}}, \code{\link[terra]{focal}}
#'   for raster I/O and preprocessing.
#'   \item \code{\link[terra]{boundaries}}, \code{\link[terra]{as.polygons}} for visualization
#'   and vectorization of segments.
#' }
#'
#' @examples
#' \dontrun{
#' library(terra)
#'
#' r <- rast("image.tif")
#'
#' seg <- fh_segmenter(
#'   r,
#'   k = 0.5,
#'   min_size = 20,
#'   eight = TRUE,
#'   scale_bands = TRUE,
#'   smooth = 3
#' )
#'
#' # Overlay boundaries on RGB
#' rgb <- stretch(r[[1:3]])
#' bnd <- boundaries(seg, directions = 8)
#' plotRGB(rgb, r=1, g=2, b=3, stretch="lin")
#' plot(bnd, add=TRUE, col="yellow", lwd=0.8)
#' }
#'
#' @export
fh_segmenter <- function(x,
                         k = 1.0,
                         min_size = 50,
                         eight = TRUE,
                         scale_bands = TRUE,
                         smooth = 0) {

  stopifnot(inherits(x, "SpatRaster"))

  # Optional smoothing (mean filter). This is safer than weighted-sum with NA omit.
  if (smooth > 0) {
    w <- matrix(1, smooth, smooth)
    x <- terra::focal(
      x,
      w = w,
      fun = mean,
      na.policy = "omit",
      fillvalue = NA
    )
  }

  v <- terra::values(x, mat = TRUE)   # ncell x nb

  # NA-aware scaling per band
  if (scale_bands) {
    for (j in seq_len(ncol(v))) {
      mu  <- mean(v[, j], na.rm = TRUE)
      sdv <- stats::sd(v[, j], na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) sdv <- 1
      v[, j] <- (v[, j] - mu) / sdv
    }
  }

  nr <- nrow(x)
  nc <- ncol(x)
  nb <- terra::nlyr(x)

  # band-major layout expected by C++: all pixels of band 1, then band 2, ...
  img <- as.numeric(v)

  lab <- fh_segmenter_cpp(
    img, nr, nc, nb,
    k = k,
    min_size = as.integer(min_size),
    eight = eight
  )

  out <- terra::rast(x, nlyr = 1)
  lab[lab == 0] <- NA_integer_
  terra::values(out) <- lab
  names(out) <- "segment_id"
  out
}
