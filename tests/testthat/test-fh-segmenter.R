test_that("fh_segmenter() basic invariants (minimal terra ops)", {
  testthat::skip_if_not_installed("terra")

  x <- sample_raster_test("memory")

  # ---- segmentation ----
  seg <- fh_segmenter(
    x,
    k = 0.8,
    min_size = 50,
    eight = TRUE,
    scale_bands = TRUE,
    smooth = 0
  )

  # ---- basic contract ----
  testthat::expect_true(inherits(seg, "SpatRaster"))
  testthat::expect_equal(terra::nlyr(seg), 1L)
  testthat::expect_true(terra::hasValues(seg))

  # geometry check (robust, avoids waldo diffs)
  testthat::expect_true(
    terra::compareGeom(seg, x, stopOnError = FALSE),
    info = "Output segmentation geometry differs from input"
  )

  # ---- "non-degenerate segmentation" checks without heavy ops ----
  # Read a small block instead of global/freq/as.matrix.
  # This confirms:
  # - values are readable
  # - at least some non-NA labels exist
  # - there is more than one label (not constant)
  v <- terra::readValues(seg, row = 1, nrows = 200, col = 1, ncols = 200, mat = FALSE)
  testthat::expect_true(length(v) > 0)

  v <- v[!is.na(v)]
  testthat::expect_gt(length(v), 0L)

  u <- unique(v)
  testthat::expect_gt(length(u), 1L)  # should not be all the same segment

  # labels should be finite numeric and mostly non-negative by convention
  testthat::expect_true(all(is.finite(u)))
  testthat::expect_gte(min(u), 0)
})

test_that("fh_segmenter() propagates NA in a small window (no as.matrix)", {
  testthat::skip_if_not_installed("terra")

  x <- sample_raster_test("file")

  # Introduce NA block in band 1
  x2 <- x
  b1 <- x2[[1]]
  b1[1:20, 1:20] <- NA
  x2[[1]] <- b1

  seg2 <- fh_segmenter(
    x2,
    k = 0.8,
    min_size = 50,
    eight = TRUE,
    scale_bands = TRUE,
    smooth = 0
  )

  # Read only the NA window from the segmentation
  vals <- terra::readValues(seg2, row = 1, nrows = 20, col = 1, ncols = 20, mat = TRUE)

  # If your intended policy is NA propagation, enforce it here
  testthat::expect_true(all(is.na(vals)))
})
