test_that("fh_meanshift_segmenter() returns a valid segmentation (core invariants)", {
  testthat::skip_if_not_installed("terra")

  x <- sample_raster_test("memory")  # helper from helper-extdata.R

  seg_ms <- fh_meanshift_segmenter(
    x,
    # preprocessing
    scale_bands = TRUE,
    smooth = 0,

    # FH stage
    fh_k = 0.7,
    fh_min_size = 30,
    eight = TRUE,

    # MeanShift stage (no PCA; ms_dim is feature dimension)
    ms_dim = 0L,         # auto: min(3, nb) is typical
    ms_ranger = 0.15,
    ms_hs = 10,
    ms_max_iter = 3L,    # keep tests fast
    ms_eps = 1e-3,

    # mode merge / cleanup
    mode_merge = 0.6,
    final_min_size = 60
  )

  # geometry + class + size invariants
  assert_seg_raster_basic(seg_ms, x)

  # label invariants (avoid expensive full scans)
  rng <- terra::global(seg_ms, range, na.rm = TRUE)
  testthat::expect_true(is.finite(rng[1, 1]))
  testthat::expect_true(is.finite(rng[1, 2]))
  testthat::expect_gte(rng[1, 1], 0)  # allow 0 if your convention uses it

  # segmentation must not be trivial (but avoid platform brittleness)
  ns <- n_segments(seg_ms)
  testthat::expect_gt(ns, 1L)
  testthat::expect_lt(ns, terra::ncell(seg_ms))
})

test_that("fh_meanshift_segmenter() is deterministic for fixed input (basic)", {
  testthat::skip_if_not_installed("terra")

  x <- sample_raster_test("memory")

  # run twice with identical parameters; should match exactly
  seg1 <- fh_meanshift_segmenter(
    x,
    scale_bands = TRUE, smooth = 0,
    fh_k = 0.7, fh_min_size = 30, eight = TRUE,
    ms_dim = 0L,
    ms_ranger = 0.15, ms_hs = 10, ms_max_iter = 3L, ms_eps = 1e-3,
    mode_merge = 0.6, final_min_size = 60
  )

  seg2 <- fh_meanshift_segmenter(
    x,
    scale_bands = TRUE, smooth = 0,
    fh_k = 0.7, fh_min_size = 30, eight = TRUE,
    ms_dim = 0L,
    ms_ranger = 0.15, ms_hs = 10, ms_max_iter = 3L, ms_eps = 1e-3,
    mode_merge = 0.6, final_min_size = 60
  )

  # Compare values, not file pointers / sources
  v1 <- terra::values(seg1, mat = FALSE)
  v2 <- terra::values(seg2, mat = FALSE)
  testthat::expect_identical(v1, v2)
})

test_that("fh_meanshift_segmenter() handles NA in the input without crashing", {
  testthat::skip_if_not_installed("terra")

  x <- sample_raster_test("memory")

  # Introduce NA into a small window (band 1 only)
  x2 <- x
  b1 <- x2[[1]]
  b1[1:20, 1:20] <- NA
  x2[[1]] <- b1

  seg_ms <- fh_meanshift_segmenter(
    x2,
    scale_bands = TRUE,
    smooth = 0,
    fh_k = 0.7,
    fh_min_size = 30,
    eight = TRUE,
    ms_dim = 0L,
    ms_ranger = 0.15,
    ms_hs = 10,
    ms_max_iter = 2L,
    ms_eps = 1e-3,
    mode_merge = 0.6,
    final_min_size = 60
  )

  assert_seg_raster_basic(seg_ms, x2)

  # If your convention is to keep invalid pixels as NA, enable this.
  # If instead you output 0 for invalid pixels, change expectation accordingly.
  # Here we only check that the region contains NA OR 0, i.e., "not positive labels".
  # Extract the same 20x20 window as a SpatRaster (not a data.frame)
  e_win <- .ext_from_rowcol(seg_ms, r0 = 1L, r1 = 20L, c0 = 1L, c1 = 20L)
  win <- terra::crop(seg_ms, e_win, snap = "out")

  vals <- terra::values(win, mat = FALSE)
  testthat::expect_true(all(is.na(vals) | vals == 0L | vals < 1L))

})
