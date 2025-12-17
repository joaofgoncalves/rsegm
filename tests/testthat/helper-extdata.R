
sample_raster_path <- function() {
  system.file("extdata", "sample_raster.tif", package = "rsegm", mustWork = TRUE)
}

sample_raster <- function() {
  terra::rast(sample_raster_path())
}

sample_raster_test <- function(materialize = c("file", "memory")) {
  materialize <- match.arg(materialize)
  x <- terra::rast(sample_raster_path())

  if (materialize == "memory") {
    x <- terra::set.values(x)
    return(x)
  }

  # file-backed materialization
  tmp <- tempfile(fileext = ".tif")
  terra::writeRaster(x, tmp, overwrite = TRUE)
  terra::rast(tmp)
}


n_segments <- function(seg) {
  fr <- terra::freq(seg, bylayer = FALSE)
  if (is.null(fr) || !nrow(fr)) return(0L)
  # ignore NA and non-positive IDs if you allow them
  v <- fr$value
  v <- v[!is.na(v) & v > 0]
  length(v)
}

assert_seg_raster_basic <- function(seg, ref) {

  expect_true(inherits(seg, "SpatRaster"))

  testthat::expect_equal(terra::nlyr(seg), 1L)

  testthat::expect_true(
    terra::compareGeom(seg, ref, stopOnError = FALSE),
    info = "seg geometry differs from reference"
  )

  testthat::expect_true(terra::hasValues(seg))
  testthat::expect_gt(terra::ncell(seg), 0L)
}
