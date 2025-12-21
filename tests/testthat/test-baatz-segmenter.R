test_that("baatz_segmenter returns a single-layer SpatRaster aligned with input", {
  skip_if_not_installed("terra")

  #x <- rast(sample_raster_path())  # returns SpatRaster

  x1 <- matrix(runif(10000, 0, 255), nrow = 100, ncol = 100)
  x2 <- matrix(runif(10000, 0, 255), nrow = 100, ncol = 100)
  x3 <- matrix(runif(10000, 0, 255), nrow = 100, ncol = 100)

  x<- c(terra::rast(x1),
        terra::rast(x2),
        terra::rast(x3))

  seg <- baatz_segmenter(x,
                         scale_param = 20,
                         color_weight = 0.5,
                         compactness = 0.5)

  expect_s4_class(seg, "SpatRaster")
  expect_equal(terra::nlyr(seg), 1L)

  # geometry alignment
  expect_equal(terra::nrow(seg), terra::nrow(x))
  expect_equal(terra::ncol(seg), terra::ncol(x))
  expect_equal(terra::crs(seg),  terra::crs(x))
  expect_identical(names(seg), "segments")

})


# test_that("baatz_segmenter is deterministic for the same inputs/parameters", {
#   skip_if_not_installed("terra")
#
#   x <- sample_raster()
#
#   seg1 <- baatz_segmenter(x, scale_param = 50, color_weight = 0.9, compactness = 0.5)
#   seg2 <- baatz_segmenter(x, scale_param = 50, color_weight = 0.9, compactness = 0.5)
#
#   v1 <- terra::values(seg1, mat = FALSE)
#   v2 <- terra::values(seg2, mat = FALSE)
#
#   expect_identical(v1, v2)
# })
#
# test_that("baatz_segmenter accepts a file path input (same result as SpatRaster)", {
#   skip_if_not_installed("terra")
#
#   x <- sample_raster()
#   f <- tempfile(fileext = ".tif")
#   terra::writeRaster(x, f, overwrite = TRUE)
#
#   seg_from_obj  <- baatz_segmenter(x)
#   seg_from_path <- baatz_segmenter(f)
#
#   expect_identical(
#     terra::values(seg_from_obj,  mat = FALSE),
#     terra::values(seg_from_path, mat = FALSE)
#   )
# })
#
# test_that("baatz_segmenter replaces NA pixels and still returns valid labels", {
#   skip_if_not_installed("terra")
#
#   x <- sample_raster()
#
#   # Introduce some NA values in the first layer only (should trigger NA handling)
#   x_na <- x
#   set.seed(1)
#   idx <- sample.int(terra::ncell(x_na), size = min(200L, terra::ncell(x_na)))
#   terra::values(x_na[[1]])[idx] <- NA
#
#   seg <- baatz_segmenter(x_na)
#
#   v <- terra::values(seg, mat = FALSE)
#   expect_false(anyNA(v))
#   expect_true(all(is.finite(v)))
#   expect_true(all(abs(v - round(v)) < 1e-8))
#   expect_gt(length(unique(v)), 1L)
# })
#
# test_that("baatz_segmenter errors on invalid input type", {
#   expect_error(baatz_segmenter(1), "Input must be a SpatRaster object or file path")
#   expect_error(baatz_segmenter(list()), "Input must be a SpatRaster object or file path")
# })
#
# test_that("baatz_segmenter errors on invalid parameter ranges (delegated to C++ is OK)", {
#   skip_if_not_installed("terra")
#
#   x <- sample_raster()
#
#   # These may be checked in C++ or may still pass; we only require they error somewhere.
#   expect_error(baatz_segmenter(x, scale_param = -1))
#   expect_error(baatz_segmenter(x, color_weight = -0.1))
#   expect_error(baatz_segmenter(x, color_weight = 1.1))
#   expect_error(baatz_segmenter(x, compactness = -0.1))
#   expect_error(baatz_segmenter(x, compactness = 1.1))
# })
#
# test_that("baatz_segmenter band_weights must match number of bands", {
#   skip_if_not_installed("terra")
#
#   x <- sample_raster()
#   nb <- terra::nlyr(x)
#
#   # correct length should work
#   bw_ok <- rep(1, nb)
#   seg <- baatz_segmenter(x, band_weights = bw_ok)
#   expect_s4_class(seg, "SpatRaster")
#
#   # incorrect length should error (in wrapper or C++)
#   bw_bad <- rep(1, nb + 1L)
#   expect_error(baatz_segmenter(x, band_weights = bw_bad))
# })
#
# test_that("baatz_segmenter writes output_file when requested", {
#   skip_if_not_installed("terra")
#
#   x <- sample_raster()
#   out <- tempfile(fileext = ".tif")
#
#   seg <- baatz_segmenter(x, output_file = out)
#
#   expect_true(file.exists(out))
#   seg2 <- terra::rast(out)
#
#   # same geometry
#   expect_equal(terra::nrow(seg2), terra::nrow(seg))
#   expect_equal(terra::ncol(seg2), terra::ncol(seg))
#   expect_equal(terra::ext(seg2),  terra::ext(seg))
#
#   # and values should match what was returned
#   expect_identical(
#     terra::values(seg,  mat = FALSE),
#     terra::values(seg2, mat = FALSE)
#   )
# })
