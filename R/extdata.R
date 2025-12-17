#' Example RGB GeoTIFF (UTM 35N)
#'
#' A small 3-band (RGB) GeoTIFF used in examples, vignettes, and tests.
#'
#' @details
#' The file is shipped with the package at:
#' `system.file("extdata", "sample_raster.tif", package = "rsegm")`.
#'
#' Read it with:
#' `terra::rast(system.file("extdata","sample_raster.tif", package="rsegm", mustWork=TRUE))`.
#'
#' @format A GeoTIFF file in `inst/extdata/`.
#' @keywords datasets
#' @name sample_raster
#'
NULL

#' Path to an example raster
#' @return A file path to an example GeoTIFF shipped with the package.
#' @export
#'
sample_raster_path <- function() {
  system.file("extdata", "sample_raster.tif", package = "rsegm", mustWork = TRUE)
}
