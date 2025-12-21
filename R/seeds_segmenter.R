
#' SEEDS superpixel segmentation
#'
#' Generates **superpixels** using a SEEDS-like hierarchical approach based on
#' fast **histogram block moves**. The method starts from a regular grid of
#' initial labels and iteratively refines boundaries by proposing block-wise
#' label changes that improve a histogram-based objective.
#'
#' This wrapper is intended to be **tile-friendly** (fast, local updates) and
#' can be used inside tiling strategies for large rasters.
#'
#' @param x A \code{\link[terra:SpatRaster-class]{SpatRaster}} with \eqn{\ge 1}
#'   layer (band).
#'
#' @param step Integer \eqn{\ge 1}. Approximate initial superpixel spacing in
#'   pixels. Larger values produce fewer (larger) superpixels. Internally this
#'   sets the size of the initial grid cells used to initialize labels.
#'
#' @param nbins Integer \eqn{\ge 2}. Number of histogram bins per band. Pixel
#'   values are quantized per band based on that band's min/max (computed
#'   ignoring \code{NA}/NaN).
#'
#' @param block_sizes Integer vector of hierarchical block sizes (coarse-to-fine),
#'   e.g. \code{c(8,4,2,1)}. At each level, candidate moves are proposed over
#'   aligned \code{bs x bs} blocks, refining boundaries progressively.
#'
#' @param iters_per_level Integer \eqn{\ge 1}. Number of refinement iterations
#'   per hierarchy level in \code{block_sizes}.
#'
#' @param boundary_samples Integer \eqn{\ge 1}. Sampling effort knob controlling
#'   how many boundary-proposal samples are evaluated per iteration.
#'   Internally the C++ routine scales the number of samples roughly with the
#'   number of labels (\code{K}) and caps it for safety.
#'
#' @param alpha Non-negative numeric. Additive smoothing constant for histogram
#'   probabilities (Laplace-like smoothing). Higher values reduce sensitivity to
#'   small counts and can stabilize results.
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
#' @param verbose Logical. If \code{TRUE}, prints progress messages from the R
#'   wrapper and allows the C++ routine to check for user interrupts during the
#'   refinement loop. If \code{FALSE}, runs quietly.
#'
#' @details
#' **Implementation notes (C++):**
#' \itemize{
#'   \item \strong{Quantization:} each band is quantized into \code{nbins} using
#'   per-band min/max computed over finite values. NaN values are quantized to
#'   \code{-1} and ignored in histogram updates.
#'   \item \strong{Initialization:} labels are initialized as a regular grid
#'   with cell size \code{step} (in pixels). Labels start as \code{0..K-1} in C++,
#'   and are returned to R as \code{1..K}.
#'   \item \strong{Moves:} at each hierarchy level, the algorithm samples boundary
#'   pixels, finds a 4-neighbour label competitor, and evaluates whether moving the
#'   aligned \code{bs x bs} block from region A to region B increases a histogram
#'   score. If beneficial, all pixels of label A inside that block are reassigned.
#'   \item \strong{Determinism:} the sampler uses a deterministic RNG in C++,
#'   yielding reproducible results given the same inputs and parameters.
#' }
#'
#' \strong{Missing values:} unlike other segmenters in \pkg{rsegm}, this SEEDS
#' implementation does not explicitly mask \code{NA} pixels. Pixels with
#' \code{NA}/NaN do not contribute to histograms, but still receive a label.
#' Consider pre-filling or masking if \code{NA} regions should be excluded.
#'
#' @return A single-layer \code{\link[terra:SpatRaster-class]{SpatRaster}} with
#'   integer superpixel IDs in a layer named \code{"seeds"}.
#'
#' @seealso \code{\link{baatz_segmenter}}, \code{\link{fh_segmenter}},
#'   \code{\link{meanshift_segmenter}}
#'
#' @examples
#' \dontrun{
#' library(terra)
#' r <- terra::rast(sample_raster_path())
#'
#' spx <- seeds_segmenter(
#'   r,
#'   step = 12,
#'   nbins = 20,
#'   block_sizes = c(8, 4, 2, 1),
#'   iters_per_level = 5,
#'   boundary_samples = 8,
#'   alpha = 1.0,
#'   verbose = TRUE
#' )
#' plot(spx)
#' }
#'
#' @export

seeds_segmenter <- function(x,
                            step = 10L,
                            nbins = 20L,
                            block_sizes = c(8L, 4L, 2L, 1L),
                            iters_per_level = 5L,
                            boundary_samples = 8L,
                            alpha = 1.0,
                            smooth=0L,
                            output_file=NULL,
                            verbose = FALSE) {

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

  # Raster -> values matrix (ncell x nb)
  v <- terra::values(x, mat = TRUE)

  if (anyNA(v)) {
    vcat("  Warning: NA values detected (handled in C++ if supported)\n")
  }

  vcat("\nRunning SEEDS superpixel segmentation...\n")
  vcat(sprintf("  Step size: %d pixels\n", step))
  vcat(sprintf("  Histogram bins: %d\n", nbins))
  vcat(sprintf("  Block sizes: %s\n", paste(block_sizes, collapse = ", ")))
  vcat(sprintf("  Iterations per level: %d\n", iters_per_level))
  vcat(sprintf("  Boundary samples: %d\n", boundary_samples))
  vcat(sprintf("  Alpha (compactness): %.3f\n", alpha))

  lab <- seeds_segmenter_cpp(
    img = v,
    nrow = nr,
    ncol = nc,
    step = as.integer(step),
    nbins = as.integer(nbins),
    block_sizes = as.integer(block_sizes),
    iters_per_level = as.integer(iters_per_level),
    boundary_samples = as.integer(boundary_samples),
    alpha = as.numeric(alpha),
    verbose = verbose
  )

  if (is.null(lab) || !is.atomic(lab)) {
    stop("seeds_superpixels_cpp() returned no labels.", call. = FALSE)
  }

  vcat("\nCreating output raster...\n")

  out <- terra::rast(x, nlyr = 1)
  terra::values(out) <- lab
  names(out) <- "seeds"

  if (isTRUE(verbose)) {
    nseg <- length(unique(lab))
    vcat("\nSegmentation complete!\n")
    vcat(sprintf("  Unique superpixels: %d\n", nseg))
    vcat(sprintf("  Average superpixel size: %.1f pixels\n",
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
