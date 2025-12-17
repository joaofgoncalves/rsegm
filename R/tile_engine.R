# ============================================================
# Utilities
# ============================================================

#' Compute maximum segment ID in a segmentation raster
#'
#' Internal utility to determine the maximum positive (non-NA) segment
#' identifier present in a segmentation raster. This is mainly used to
#' offset segment IDs safely when merging tiled segmentations.
#'
#' @param seg A `SpatRaster` containing integer segment labels.
#'
#' @return An integer scalar giving the maximum segment ID found,
#'   or `0L` if no valid segments are present.
#'
#' @keywords internal
#' @noRd
.tile_max_id <- function(seg) {
  fr <- terra::freq(seg, bylayer = FALSE)
  if (is.null(fr) || !nrow(fr)) return(0L)
  v <- fr$value
  v <- v[!is.na(v) & v > 0]
  if (!length(v)) return(0L)
  as.integer(max(v))
}

#' Test whether a raster is readable
#'
#' Internal helper that checks whether a `SpatRaster` can be read without
#' error by attempting to access a small block of values. This is used
#' defensively in tiled and streaming workflows to detect invalid or
#' corrupted rasters early.
#'
#' @param r A `SpatRaster` object.
#'
#' @return Logical scalar. `TRUE` if the raster can be read successfully,
#'   `FALSE` otherwise.
#'
#' @keywords internal
#' @noRd
.is_readable_raster <- function(r) {
  tryCatch({
    nr <- terra::nrow(r)
    if (is.na(nr) || nr < 1) return(FALSE)
    v <- terra::readValues(r, row = 1, nrows = min(8L, nr), mat = FALSE)
    !is.null(v) && length(v) > 0
  }, error = function(e) FALSE)
}


#' Default GDAL creation options
#'
#' Internal helper returning a standard set of GDAL creation options used
#' when writing raster outputs to disk. The defaults favor compression,
#' tiling, and safe BigTIFF handling for large segmentation products.
#'
#' @return A character vector of GDAL creation options.
#'
#' @keywords internal
#' @noRd
.default_gdal_opts <- function() c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER")

#' Convert row/column window to spatial extent
#'
#' Internal utility that converts a raster row/column window into a spatial
#' extent (`SpatExtent`). The computation is resolution-aware and includes
#' grid snapping to avoid rounding artifacts along tile borders. Designed
#' to be safe with terra >= 1.8-60.
#'
#' @param x A reference `SpatRaster`.
#' @param r0,r1 Integer row indices (start and end, inclusive).
#' @param c0,c1 Integer column indices (start and end, inclusive).
#'
#' @return A `SpatExtent` corresponding to the specified window.
#'
#' @keywords internal
#' @noRd
# row/col window -> spatial extent (terra 1.8-60 safe)
.ext_from_rowcol <- function(x, r0, r1, c0, c1) {
  e  <- terra::ext(x)
  rs <- terra::res(x)
  xr <- rs[1]
  yr <- rs[2]

  r0 <- as.integer(r0)
  r1 <- as.integer(r1)

  c0 <- as.integer(c0)
  c1 <- as.integer(c1)

  xmin_win <- terra::xmin(e) + (c0 - 1L) * xr
  xmax_win <- terra::xmin(e) + (c1)      * xr
  ymax_win <- terra::ymax(e) - (r0 - 1L) * yr
  ymin_win <- terra::ymax(e) - (r1)      * yr

  # snap to grid to avoid border rounding empties
  xmin_win <- terra::xmin(e) + round((xmin_win - terra::xmin(e)) / xr) * xr
  xmax_win <- terra::xmin(e) + round((xmax_win - terra::xmin(e)) / xr) * xr
  ymin_win <- terra::ymax(e) - round((terra::ymax(e) - ymin_win) / yr) * yr
  ymax_win <- terra::ymax(e) - round((terra::ymax(e) - ymax_win) / yr) * yr

  terra::ext(c(xmin_win, xmax_win, ymin_win, ymax_win))
}

# ============================================================
# 1) Deterministic tiling to disk
# ============================================================

#' Create deterministic raster tiles on disk
#'
#' Splits a large `SpatRaster` into fixed-size tiles with optional overlap
#' (buffer) and writes each tile to disk. Tiles are generated deterministically
#' in row-column order and accompanied by metadata describing their spatial
#' layout and buffered extents.
#'
#' This function is designed for scalable, disk-backed segmentation workflows
#' where large rasters must be processed tile-by-tile while preserving spatial
#' consistency at tile borders.
#'
#' @param x A `SpatRaster` to be tiled.
#' @param tile_size Integer. Size (in pixels) of each tile along rows and
#'   columns (excluding buffer).
#' @param buffer Integer. Number of pixels added as overlap on all sides of
#'   each tile.
#' @param prefix Character. Filename prefix for generated tile files.
#' @param dir Character. Output directory where tiles will be written.
#' @param overwrite Logical. Whether to overwrite existing tile files.
#' @param verbose Logical. If `TRUE`, display progress information.
#'
#' @return A `data.frame` with one row per tile and the following columns:
#'   \describe{
#'     \item{tile_id}{Sequential tile identifier.}
#'     \item{file}{Path to the tile raster on disk.}
#'     \item{r0_in, r1_in}{Row indices of the inner (unbuffered) tile.}
#'     \item{c0_in, c1_in}{Column indices of the inner (unbuffered) tile.}
#'     \item{r0_buf, r1_buf}{Row indices including buffer.}
#'     \item{c0_buf, c1_buf}{Column indices including buffer.}
#'     \item{buffer}{Buffer size (in pixels).}
#'   }
#'
#' @details
#' Tiles are cropped using spatial extents derived from row/column indices,
#' snapped to the raster grid to avoid rounding artifacts. Output rasters are
#' written using tiled GeoTIFF with LZW compression and safe BigTIFF handling.
#'
#' Empty tiles (with no values) are skipped. An error is raised if no tiles are
#' produced.
#'
#' @seealso
#' \code{\link[terra]{crop}}, \code{\link[terra]{writeRaster}}
#'
#' @export
#'
make_tiles_disk <- function(x, tile_size=2048, buffer=64, prefix="tile",
                            dir=tempdir(), overwrite=TRUE, verbose=TRUE) {

  stopifnot(inherits(x, "SpatRaster"))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  nr <- terra::nrow(x)
  nc <- terra::ncol(x)

  tile_size <- as.integer(tile_size)
  buffer <- as.integer(buffer)

  r_starts <- seq.int(1L, nr, by = tile_size)
  c_starts <- seq.int(1L, nc, by = tile_size)

  t_len <- length(r_starts) * length(c_starts)

  out <- list()
  k <- 0L

  if(verbose){
    message("\n## [1/7] Creating tiles:\n")
    pb <- utils::txtProgressBar(min = 1, max = t_len, style = 3)
  }


  for (r0 in r_starts) {
    r1_in <- min(nr, r0 + tile_size - 1L)
    for (c0 in c_starts) {
      c1_in <- min(nc, c0 + tile_size - 1L)

      r0_buf <- max(1L, r0 - buffer)
      c0_buf <- max(1L, c0 - buffer)
      r1_buf <- min(nr, r1_in + buffer)
      c1_buf <- min(nc, c1_in + buffer)

      e <- .ext_from_rowcol(x, r0_buf, r1_buf, c0_buf, c1_buf)
      tile <- terra::crop(x, e, snap="out")

      if (!terra::hasValues(tile) || terra::ncell(tile) == 0) next

      k <- k + 1L
      f <- file.path(dir, sprintf("%s_%05d.tif", prefix, k))
      if (file.exists(f) && overwrite) file.remove(f)

      terra::writeRaster(tile, f, overwrite=TRUE, gdal=.default_gdal_opts())

      out[[k]] <- data.frame(
        tile_id=k, file=f,
        r0_in=r0, r1_in=r1_in, c0_in=c0, c1_in=c1_in,
        r0_buf=r0_buf, r1_buf=r1_buf, c0_buf=c0_buf, c1_buf=c1_buf,
        buffer=buffer, stringsAsFactors=FALSE
      )

      if(verbose) utils::setTxtProgressBar(pb, k)

    }
  }

  if (!length(out)) stop("No tiles were produced (all tiles empty?)")
  do.call(rbind, out)
}

# ============================================================
# 2) Seam mask / seam adjacency
# ============================================================

#' Build a logical mask of seam pixels for a buffered tile
#'
#' Internal helper that constructs a logical matrix identifying the seam zone
#' (overlap/buffer area) around the *inner* portion of a buffered tile. The seam
#' mask is typically used to restrict adjacency/pair extraction to pixels near
#' tile borders when reconciling segment IDs across tiles.
#'
#' @param seg_tile A `SpatRaster` tile segmentation (usually a buffered tile)
#'   with integer segment labels.
#' @param tile_meta_row A single-row `data.frame` (or list-like object)
#'   containing the tiling metadata for the tile, as returned by
#'   [make_tiles_disk()]. Must include `r0_in`, `r1_in`, `c0_in`, `c1_in`,
#'   `r0_buf`, `r1_buf`, `c0_buf`, `c1_buf`, and `buffer`.
#' @param seam_width Integer. Width (in pixels) of the seam zone to mark. If
#'   `NULL`, defaults to `tile_meta_row$buffer`.
#'
#' @return A logical matrix with dimensions `nrow(seg_tile)` by `ncol(seg_tile)`,
#'   where `TRUE` indicates pixels considered part of the seam zone.
#'
#' @details
#' If `seam_width <= 0`, the function returns an all-`TRUE` mask. The effective
#' seam width is clipped to avoid consuming the inner tile entirely (useful for
#' small tiles or large buffers).
#'
#' @keywords internal
#' @noRd
#'
make_tile_seam_mask <- function(seg_tile, tile_meta_row, seam_width = NULL) {

  stopifnot(inherits(seg_tile, "SpatRaster"))

  if (is.null(seam_width)) seam_width <- tile_meta_row$buffer
  b <- as.integer(seam_width)

  nr <- terra::nrow(seg_tile)
  nc <- terra::ncol(seg_tile)

  if (b <= 0L) return(matrix(TRUE, nr, nc))

  top_off   <- (tile_meta_row$r0_in - tile_meta_row$r0_buf) + 1L
  left_off  <- (tile_meta_row$c0_in - tile_meta_row$c0_buf) + 1L
  bot_off   <- (tile_meta_row$r1_in - tile_meta_row$r0_buf) + 1L
  right_off <- (tile_meta_row$c1_in - tile_meta_row$c0_buf) + 1L

  inner_h <- bot_off - top_off + 1L
  inner_w <- right_off - left_off + 1L
  b <- min(b, floor(min(inner_h, inner_w)/2) - 1L)
  if (b < 1L) return(matrix(TRUE, nr, nc))

  m <- matrix(FALSE, nr, nc)

  r_top1 <- top_off
  r_top2 <- min(nr, top_off + b - 1L)
  r_bot1 <- max(1L, bot_off - b + 1L)
  r_bot2 <- bot_off

  c_left1  <- left_off
  c_left2  <- min(nc, left_off + b - 1L)
  c_right1 <- max(1L, right_off - b + 1L)
  c_right2 <- right_off

  m[r_top1:r_top2, ] <- TRUE
  m[r_bot1:r_bot2, ] <- TRUE
  m[, c_left1:c_left2] <- TRUE
  m[, c_right1:c_right2] <- TRUE

  m
}

#' Extract unique adjacency pairs of segment IDs along seams
#'
#' Internal helper that computes unique neighboring segment label pairs
#' (adjacencies) from a segmentation raster. Optionally restricts computation
#' to a seam zone (e.g., the tile overlap area) via a logical mask.
#'
#' Adjacent pairs are returned as undirected edges: each pair is sorted
#' `(min, max)` and duplicates are removed.
#'
#' @param seg A `SpatRaster` with integer segment labels.
#' @param seam_mask Optional logical matrix with dimensions matching `seg`.
#'   If provided, only adjacencies where at least one of the two neighboring
#'   pixels lies in the masked area are retained. This is typically produced by
#'   [make_tile_seam_mask()].
#' @param directions Integer. Neighborhood connectivity: `4` (rook) or `8`
#'   (queen).
#'
#' @return An integer matrix with two columns, where each row is a unique pair
#'   of adjacent segment IDs. Returns a `0 x 2` integer matrix if no valid pairs
#'   are found.
#'
#' @details
#' Segment IDs that are `NA` are ignored. Self-pairs (same ID on both sides) are
#' discarded. When `seam_mask` is provided, the raster may be cropped internally
#' to the bounding region of the seam to reduce work.
#'
#' @keywords internal
#' @noRd
#'
extract_seam_pairs <- function(seg, seam_mask = NULL, directions = 8) {

  stopifnot(inherits(seg, "SpatRaster"))
  stopifnot(directions %in% c(4, 8))

  if (!is.null(seam_mask)) {
    if (!is.matrix(seam_mask)) stop("seam_mask must be a logical matrix.")

    nr0 <- terra::nrow(seg)
    nc0 <- terra::ncol(seg)

    if (any(dim(seam_mask) != c(nr0, nc0))) stop("seam_mask dims != seg dims")

    idx <- which(seam_mask, arr.ind = TRUE)
    if (!nrow(idx)) return(matrix(integer(0), ncol = 2))

    rmin <- max(1L, min(idx[, 1]) - 1L)
    rmax <- min(nr0, max(idx[, 1]) + 1L)
    cmin <- max(1L, min(idx[, 2]) - 1L)
    cmax <- min(nc0, max(idx[, 2]) + 1L)

    e <- .ext_from_rowcol(seg, rmin, rmax, cmin, cmax)
    seg <- terra::crop(seg, e, snap = "out")
    seam_mask <- seam_mask[rmin:rmax, cmin:cmax, drop = FALSE]
  }

  m <- terra::as.matrix(seg, wide = TRUE)
  nr <- nrow(m)
  nc <- ncol(m)

  collect_pairs <- function(a, b, keep_mask = NULL) {
    ok <- !is.na(a) & !is.na(b) & (a != b)
    if (!is.null(keep_mask)) ok <- ok & keep_mask
    if (!any(ok)) return(NULL)
    cbind(a[ok], b[ok])
  }

  pairs <- NULL

  a <- m[, 1:(nc-1), drop = FALSE]
  b <- m[, 2:nc,     drop = FALSE]
  keep <- NULL

  if (!is.null(seam_mask)){
    keep <- seam_mask[, 1:(nc-1), drop=FALSE] | seam_mask[, 2:nc, drop=FALSE]
  }

  p <- collect_pairs(a, b, keep); if (!is.null(p)) pairs <- rbind(pairs, p)

  a <- m[1:(nr-1), , drop = FALSE]
  b <- m[2:nr,     , drop = FALSE]
  keep <- NULL

  if (!is.null(seam_mask)){
    keep <- seam_mask[1:(nr-1), , drop=FALSE] | seam_mask[2:nr, , drop=FALSE]
  }

  p <- collect_pairs(a, b, keep); if (!is.null(p)) pairs <- rbind(pairs, p)

  if (directions == 8) {

    a <- m[1:(nr-1), 1:(nc-1), drop = FALSE]
    b <- m[2:nr,     2:nc,     drop = FALSE]
    keep <- NULL

    if (!is.null(seam_mask)){
      keep <- seam_mask[1:(nr-1), 1:(nc-1), drop=FALSE] | seam_mask[2:nr, 2:nc, drop=FALSE]
    }
    p <- collect_pairs(a, b, keep); if (!is.null(p)) pairs <- rbind(pairs, p)

    a <- m[1:(nr-1), 2:nc,     drop = FALSE]
    b <- m[2:nr,     1:(nc-1), drop = FALSE]
    keep <- NULL

    if (!is.null(seam_mask)){
      keep <- seam_mask[1:(nr-1), 2:nc, drop=FALSE] | seam_mask[2:nr, 1:(nc-1), drop=FALSE]
    }
    p <- collect_pairs(a, b, keep); if (!is.null(p)) pairs <- rbind(pairs, p)
  }

  if (is.null(pairs) || !nrow(pairs)) return(matrix(integer(0), ncol = 2))

  pairs <- t(apply(pairs, 1, sort))
  pairs <- unique(pairs)
  pairs[pairs[,1] != pairs[,2], , drop = FALSE]
}

# ============================================================
# 3) Means from tiles (no VRT reads)
# ============================================================

#' Compute per-segment band means from tiled segmentations
#'
#' Computes mean spectral values for a set of segment IDs by iterating over
#' segmentation tiles stored on disk and streaming values from the original
#' image raster. This function avoids virtual raster (VRT) reads and is designed
#' for memory-safe processing of large images.
#'
#' @param img A `SpatRaster` containing the original multi-band image data.
#' @param seg_tile_files Character vector of file paths to segmentation tiles
#'   (e.g., produced by [make_tiles_disk()]).
#' @param ids Integer vector of segment IDs for which means should be computed.
#' @param block_nrows Integer. Number of raster rows to process per read block
#'   when streaming values from disk.
#'
#' @return A numeric matrix with `length(ids)` rows and `nlyr(img)` columns.
#'   Rows correspond to segment IDs (in the order of `ids`), and columns
#'   correspond to image bands. Segments with no contributing pixels return
#'   `NA` values.
#'
#' @details
#' For each segmentation tile, the corresponding spatial subset of `img` is
#' cropped and read in row blocks. If grid geometry differs between the image
#' and a tile, the image subset is resampled to the tile grid before streaming.
#'
#' Accumulation is performed as sums and pixel counts per segment, followed by
#' normalization to means. Segment labels with `NA` values are ignored.
#'
#' This function is typically used after tiled segmentation and ID reconciliation
#' to compute region statistics without loading full rasters into memory.
#'
#' @seealso
#' \code{\link[terra]{readValues}}, \code{\link[terra]{crop}},
#' \code{\link[terra]{resample}}
#'
#' @export
#'
segment_means_from_tiles <- function(img, seg_tile_files, ids, block_nrows=512) {

  stopifnot(inherits(img, "SpatRaster"))
  ids <- as.integer(ids)
  nb  <- terra::nlyr(img)

  sum <- matrix(0, nrow=length(ids), ncol=nb, dimnames=list(as.character(ids), NULL))
  cnt <- integer(length(ids))

  for (f in seg_tile_files) {
    segt <- terra::rast(f)
    if (!terra::hasValues(segt) || terra::ncell(segt) == 0) next

    if (!tryCatch({ terra::readStart(segt); TRUE }, error=function(e) FALSE)) next

    imgt <- terra::crop(img, terra::ext(segt), snap="out")

    # ensure readable / concrete if lazy
    if (!terra::hasValues(imgt) || terra::ncell(imgt) == 0) {
      terra::readStop(segt)
      next
    }

    # align grid if needed (materialize resample)
    if (terra::nrow(imgt) != terra::nrow(segt) ||
        terra::ncol(imgt) != terra::ncol(segt) ||
        any(terra::res(imgt) != terra::res(segt))) {
      tmp_rs <- tempfile(fileext = ".tif")
      rrs <- terra::resample(imgt, segt, method="bilinear")
      terra::writeRaster(rrs, tmp_rs, overwrite=TRUE, gdal=.default_gdal_opts())
      imgt <- terra::rast(tmp_rs)
    }

    if (!tryCatch({ terra::readStart(imgt); TRUE }, error=function(e) FALSE)) {
      terra::readStop(segt)
      next
    }

    nr <- terra::nrow(segt)
    for (r0 in seq.int(1L, nr, by=as.integer(block_nrows))) {
      nr_here <- min(as.integer(block_nrows), nr - r0 + 1L)

      lab <- terra::readValues(segt, row=r0, nrows=nr_here, mat=FALSE)
      if (is.null(lab) || !length(lab)) next

      dat <- terra::readValues(imgt, row=r0, nrows=nr_here, mat=TRUE)
      if (is.null(dat) || !nrow(dat)) next

      ok <- !is.na(lab)
      if (!any(ok)) next

      lab_ok <- lab[ok]
      dat_ok <- dat[ok, , drop=FALSE]

      idx <- match(lab_ok, ids)
      keep <- !is.na(idx)
      if (!any(keep)) next

      idx <- idx[keep]
      dat_ok <- dat_ok[keep, , drop=FALSE]

      rs <- rowsum(dat_ok, idx, reorder=FALSE)
      grp <- as.integer(rownames(rs))
      sum[grp, ] <- sum[grp, , drop=FALSE] + rs

      cnt <- cnt + tabulate(idx, nbins=length(ids))
    }

    try(terra::readStop(segt), silent=TRUE)
    try(terra::readStop(imgt), silent=TRUE)
  }

  means <- sum
  nz <- cnt > 0L
  means[nz, ] <- means[nz, , drop=FALSE] / cnt[nz]
  means[!nz, ] <- NA_real_
  means
}

# ============================================================
# 4) merge_by_adjacency
# ============================================================

#' Merge adjacent segments by spectral similarity
#'
#' Given a set of adjacent segment ID pairs and per-segment mean feature vectors,
#' compute a mapping that merges (clusters) segments whose mean vectors are
#' sufficiently similar. Similarity is evaluated only between segments that are
#' adjacent (share an edge/corner), and merges are propagated transitively via
#' a union-find (disjoint-set) structure.
#'
#' This function is typically used after extracting seam or boundary adjacencies
#' (e.g., from tiled segmentations) to reconcile segment IDs across tiles and
#' reduce seam artifacts by merging regions with near-identical spectral means.
#'
#' @param adj_pairs Integer matrix with two columns. Each row is an undirected
#'   adjacency between segment IDs (e.g., output of `extract_seam_pairs()`).
#'   Pairs may contain duplicates; they do not need to be sorted.
#' @param means Numeric matrix of per-segment mean features (e.g., band means).
#'   Row names must be segment IDs as character strings (e.g., `"123"`), and
#'   columns are feature dimensions (bands or derived features).
#' @param thr Numeric scalar. Merge threshold applied to the Euclidean distance
#'   between mean feature vectors. Adjacent segments with distance `< thr` are
#'   merged into the same component.
#'
#' @return An integer vector mapping segment IDs to representative IDs. The
#'   vector is named by the original segment IDs (as character) and its values
#'   are the representative IDs (as integer). Only IDs present in `adj_pairs`
#'   are included. If `adj_pairs` is empty, returns an empty named integer
#'   vector.
#'
#' @details
#' **Algorithm.** For each adjacency `(a, b)`, the function extracts feature
#' vectors `means[a, ]` and `means[b, ]` and computes their Euclidean distance
#' `d = sqrt(sum((da - db)^2))`. If `d < thr` and both vectors contain no `NA`,
#' the segments are unioned in a disjoint-set (union-find) structure. After all
#' pairs are processed, each segment is assigned to the representative of its
#' connected component (with path compression to accelerate repeated finds).
#'
#' **Transitivity.** Because union-find forms connected components, merges are
#' transitive: if `a` merges with `b` and `b` merges with `c`, then `a`, `b`,
#' and `c` will share the same representative even if `(a, c)` was never an
#' explicit adjacency.
#'
#' **Representatives.** Representatives are chosen by union order (the root in
#' the union-find structure), not necessarily the smallest ID. If you require
#' canonical representatives (e.g., minimum ID per component), post-process the
#' mapping accordingly.
#'
#' **Robustness.** Pairs whose feature vectors contain `NA` are skipped. Only
#' IDs appearing in `adj_pairs` are considered; if `means` contains additional
#' segments not present in `adj_pairs`, they will not appear in the output map.
#'
#' @examples
#' \dontrun{
#'
#' # Example adjacency pairs (undirected)
#' adj <- matrix(c(1,2, 2,3, 10,11), ncol = 2, byrow = TRUE)
#' means <- rbind(
#'   "1"  = c(0.10, 0.20),
#'   "2"  = c(0.11, 0.19),
#'   "3"  = c(0.12, 0.18),
#'   "10" = c(0.80, 0.75),
#'   "11" = c(0.81, 0.74)
#' )
#' merge_by_adjacency(adj, means, thr = 0.05)
#' }
#'
#' @export
#'
merge_by_adjacency <- function(adj_pairs, means, thr) {
  # adj_pairs: integer matrix (n x 2)
  # means: matrix with rownames = segment ids as characters

  if (is.null(adj_pairs) || !nrow(adj_pairs)) return(stats::setNames(integer(0), character(0)))

  ids <- sort(unique(as.integer(c(adj_pairs))))
  parent <- seq_along(ids)

  find <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }

  unite <- function(i, j) {
    ri <- find(i)
    rj <- find(j)
    if (ri != rj) parent[rj] <<- ri
  }

  for (k in seq_len(nrow(adj_pairs))) {
    a <- adj_pairs[k, 1]
    b <- adj_pairs[k, 2]

    da <- means[as.character(a), , drop=TRUE]
    db <- means[as.character(b), , drop=TRUE]

    if (anyNA(da) || anyNA(db)) next

    d <- sqrt(sum((da - db)^2))
    if (is.finite(d) && d < thr) {
      ia <- match(a, ids)
      ib <- match(b, ids)
      if (!is.na(ia) && !is.na(ib)) unite(ia, ib)
    }
  }

  rep <- vapply(seq_along(ids), find, integer(1))
  map <- ids[rep]
  stats::setNames(as.integer(map), as.character(ids))
}

# ============================================================
# 5) Apply merge map + global relabel (streaming-safe)
# ============================================================

#' Apply a segment merge map to a segmentation raster (streaming-safe)
#'
#' Applies a mapping from original segment IDs to merged (representative) IDs
#' and writes the updated segmentation to disk using block-wise streaming.
#' This is intended for large, file-backed rasters and avoids reading the full
#' segmentation into memory.
#'
#' @param seg A `SpatRaster` containing integer segment labels.
#' @param map Named integer vector defining the relabeling/merge map. Names are
#'   original segment IDs (as character strings) and values are the target IDs
#'   (typically representatives produced by `merge_by_adjacency()`).
#' @param out_file Character. Output file path (GeoTIFF recommended).
#' @param block_nrows Integer. Number of rows per streaming block when reading
#'   and writing.
#'
#' @return A file-backed `SpatRaster` pointing to `out_file`.
#'
#' @details
#' The function streams through `seg` in row blocks, replaces values that match
#' `names(map)`, and writes the result to `out_file`. Values not present in the
#' map are left unchanged. `NA` values are preserved.
#'
#' **terra I/O note.** In terra, `writeStart()` returns a list; this function
#' intentionally ignores that return value and continues writing to the
#' `SpatRaster` object used in `writeStart()`. This ensures compatibility with
#' terra versions where writing to the returned list is incorrect.
#'
#' The output is written with integer datatype (`INT4S`) and default GDAL
#' creation options used by the package.
#'
#' @seealso
#' `merge_by_adjacency()`, `relabel_segments_global()`,
#' \code{\link[terra]{readValues}}, \code{\link[terra]{writeValues}}
#'
#' @export
#'
apply_merge_map <- function(seg, map, out_file, block_nrows = 1024L) {
  stopifnot(inherits(seg, "SpatRaster"))
  if (!terra::hasValues(seg)) stop("apply_merge_map(): seg has no values.")

  nr <- terra::nrow(seg)
  block_nrows <- as.integer(block_nrows)
  if (is.na(block_nrows) || block_nrows < 1L) block_nrows <- 1024L

  keys <- as.integer(names(map))
  vals <- as.integer(unname(map))

  # open seg for reading
  terra::readStart(seg)
  on.exit(try(terra::readStop(seg), silent = TRUE), add = TRUE)

  out <- terra::rast(seg)

  # IMPORTANT: writeStart returns a LIST in terra, but you IGNORE it.
  invisible(terra::writeStart(
    out, out_file, overwrite = TRUE,
    wopt = list(datatype = "INT4S"),
    gdal = .default_gdal_opts()
  ))
  on.exit(try(terra::writeStop(out), silent = TRUE), add = TRUE)

  for (r0 in seq.int(1L, nr, by = block_nrows)) {
    nr_here <- min(block_nrows, nr - r0 + 1L)

    v <- terra::readValues(seg, row = r0, nrows = nr_here, mat = FALSE)
    if (is.null(v) || !length(v)) next

    hit <- match(v, keys)
    ok  <- !is.na(hit)
    if (any(ok)) v[ok] <- vals[hit[ok]]

    # writeValues must receive the SpatRaster 'out', NOT the list from writeStart
    terra::writeValues(out, v, start = r0, nrows = nr_here)
  }

  terra::writeStop(out)
  terra::rast(out_file)
}

#' Relabel segments to a global consecutive ID scheme (streaming-safe)
#'
#' Reassigns all non-`NA` segment IDs in a segmentation raster to a consecutive
#' integer sequence `1..K` and writes the result to disk using streaming I/O.
#' This is useful after merging/reconciliation steps to produce compact,
#' stable IDs for downstream feature extraction and polygonization.
#'
#' @param seg A `SpatRaster` containing integer segment labels.
#' @param out_file Character. Output file path (GeoTIFF recommended).
#' @param block_nrows Integer. Number of rows per streaming block when reading
#'   and writing.
#'
#' @return A file-backed `SpatRaster` pointing to `out_file`.
#'
#' @details
#' The function first ensures `seg` is file-backed (written to a temporary
#' GeoTIFF) to support robust streaming reads. It then computes the set of
#' unique segment IDs using `terra::freq(..., bylayer = FALSE)` and maps them in
#' sorted order to `1..K`. `NA` values are ignored and preserved.
#'
#' The relabeling is applied block-wise using `readValues()` / `writeValues()`,
#' avoiding full in-memory materialization. Output is written as `INT4S` with
#' the package's default GDAL creation options.
#'
#' **Determinism.** Because unique IDs are sorted before mapping, the resulting
#' labels are deterministic for a given raster (subject to identical numeric
#' values in `seg`).
#'
#' @seealso
#' `apply_merge_map()`, \code{\link[terra]{freq}},
#' \code{\link[terra]{readValues}}, \code{\link[terra]{writeValues}}
#'
#' @export
#'
relabel_segments_global <- function(seg, out_file, block_nrows = 1024L) {
  # Ensure file-backed for stable streaming
  tmp <- tempfile(fileext = ".tif")
  terra::writeRaster(seg, tmp, overwrite = TRUE,
                     wopt = list(datatype = "INT4S"),
                     gdal = .default_gdal_opts())
  seg <- terra::rast(tmp)

  # terra 1.8-60 compatible: bylayer=FALSE returns (value,count)
  fr <- terra::freq(seg, bylayer = FALSE)
  if (is.null(fr) || !nrow(fr)) stop("relabel_segments_global(): freq() returned empty table.")
  u <- sort(as.integer(fr$value[!is.na(fr$value)]))
  if (!length(u)) stop("relabel_segments_global(): no non-NA segment IDs found.")

  keys <- u
  vals <- seq_along(u)

  nr <- terra::nrow(seg)
  block_nrows <- as.integer(block_nrows)
  if (is.na(block_nrows) || block_nrows < 1L) block_nrows <- 1024L

  terra::readStart(seg)
  on.exit(try(terra::readStop(seg), silent = TRUE), add = TRUE)

  out <- terra::rast(seg)

  # NOTE: writeStart returns a list in terra; ignore it and keep writing to 'out'
  invisible(terra::writeStart(
    out, out_file, overwrite = TRUE,
    wopt = list(datatype = "INT4S"),
    gdal = .default_gdal_opts()
  ))
  on.exit(try(terra::writeStop(out), silent = TRUE), add = TRUE)

  for (r0 in seq.int(1L, nr, by = block_nrows)) {
    nr_here <- min(block_nrows, nr - r0 + 1L)
    v <- terra::readValues(seg, row = r0, nrows = nr_here, mat = FALSE)
    if (is.null(v) || !length(v)) next

    hit <- match(v, keys)
    ok  <- !is.na(hit)
    if (any(ok)) v[ok] <- vals[hit[ok]]

    terra::writeValues(out, v, start = r0, nrows = nr_here)
  }

  terra::writeStop(out)
  terra::rast(out_file)
}


# ============================================================
# 6) Segment each tile to disk
# ============================================================

#' Segment each raster tile to disk with globally unique IDs
#'
#' Applies a user-provided segmentation function to each tile described in a
#' tiling metadata table and writes the resulting 1-layer integer segmentation
#' to disk. Segment IDs are made globally unique across tiles by applying a
#' running offset to positive labels.
#'
#' This function is designed to be used as a step in the tiled segmentation
#' engine, but can also be used standalone when users want explicit control
#' over tiling/segmentation and subsequent merging.
#'
#' @param tile_meta A `data.frame` describing input tiles on disk. Must contain
#'   a `file` column with paths to tile rasters (e.g., output of
#'   [make_tiles_disk()]).
#' @param segment_fun A function implementing segmentation on a single tile.
#'   It must accept a `SpatRaster` as first argument and return a 1-layer
#'   `SpatRaster` of integer segment labels (or something coercible to that).
#' @param seg_args A named list of additional arguments passed to `segment_fun`
#'   via `do.call()`.
#' @param out_dir Character. Directory where segmented tile rasters will be
#'   written.
#' @param out_prefix Character. Filename prefix for segmented tile outputs.
#' @param overwrite Logical. Whether to overwrite existing output files.
#' @param cleanup_tiles Logical. If `TRUE`, delete the input tile rasters after
#'   successful segmentation.
#' @param force_positive Logical. If `TRUE`, assumes the segmenter may emit
#'   non-positive labels and enforces semantics where only labels `> 0` are
#'   offset to ensure global uniqueness (leaving `0` and `NA` unchanged).
#' @param verbose Logical. If `TRUE`, display progress information.
#'
#' @return A character vector of file paths to segmented tile rasters (one per
#'   row in `tile_meta`), in processing order.
#'
#' @details
#' **Failure handling.** If segmentation fails on a tile (error or no values),
#' an all-`NA` segmentation tile is written to disk so downstream steps can
#' proceed deterministically.
#'
#' **Global uniqueness.** After each tile is segmented, the maximum positive ID
#' is computed and used to update a global offset. Positive labels are shifted
#' by the current offset to prevent ID collisions across tiles. The output is
#' written as `INT4S` (32-bit signed integer); an error is raised if the offset
#' would overflow this range.
#'
#' @seealso
#' [make_tiles_disk()], `.tile_max_id()`, \code{\link[terra]{writeRaster}}
#'
#' @export
#'
segment_tiles <- function(tile_meta, segment_fun, seg_args,
                          out_dir,
                          out_prefix = "seg",
                          overwrite = TRUE,
                          cleanup_tiles = FALSE,
                          force_positive = TRUE,
                          verbose=TRUE) {

  stopifnot(is.data.frame(tile_meta), "file" %in% names(tile_meta))
  stopifnot(!missing(out_dir))
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  out_files <- character(nrow(tile_meta))

  if(verbose){
    message("\n## [2/7] Performing segmentation on each tile:\n")
    pb <- utils::txtProgressBar(min = 1, max = nrow(tile_meta), style = 3)
  }

  global_offset <- 0L

  for (i in seq_len(nrow(tile_meta))) {
    tile_path <- tile_meta$file[i]
    r <- terra::rast(tile_path)

    out_file <- file.path(out_dir, sprintf("%s_%05d.tif", out_prefix, i))
    tmp_file <- file.path(out_dir, sprintf("%s_%05d_tmp.tif", out_prefix, i))

    ok <- TRUE
    seg <- tryCatch(
      do.call(segment_fun, c(list(r), seg_args)),
      error = function(e) {
        ok <<- FALSE
        message(sprintf("Segmentation failed on tile %d (%s): %s", i, tile_path, e$message))
        NULL
      }
    )

    if (!ok || is.null(seg) || !terra::hasValues(seg)) {
      seg <- terra::rast(r, nlyr = 1)
      terra::values(seg) <- NA_integer_
      names(seg) <- "segment_id"
    }

    # Ensure 1-layer integer labels
    if (terra::nlyr(seg) != 1) seg <- seg[[1]]
    names(seg) <- "segment_id"

    # Optional: if your segmenter can output 0/-ve labels, normalize
    if (force_positive) {
      # Keep NA as NA, keep 0 as 0, shift only positives later
      # (If you want 0 to become NA, do it here explicitly.)
      # seg[seg == 0] <- NA
      NULL
    }

    # ---- CRITICAL FIX ----
    # Make IDs globally unique by offsetting positive labels
    max_local <- .tile_max_id(seg)

    if (max_local > 0L) {
      # Apply offset only to >0 labels (preserve 0/NA semantics)
      seg <- terra::ifel(!is.na(seg) & seg > 0, seg + global_offset, seg)

      # Update offset for next tile
      # Guard against INT32 overflow (GeoTIFF INT4S)
      if ((global_offset + max_local) > 2147483647L) {
        stop("Segment IDs exceed INT32 range. Use a 64-bit datatype strategy.")
      }
      global_offset <- global_offset + max_local
    }

    terra::writeRaster(
      seg, tmp_file, overwrite = TRUE,
      wopt = list(datatype = "INT4S"),
      gdal = .default_gdal_opts()
    )

    if (file.exists(out_file) && overwrite) file.remove(out_file)
    file.rename(tmp_file, out_file)

    out_files[i] <- out_file

    if (cleanup_tiles) try(file.remove(tile_path), silent = TRUE)
    gc()
    if(verbose) utils::setTxtProgressBar(pb, i)
  }

  out_files
}



# ============================================================
# 7) Main tiled engine (clean, no duplicated blocks)
# ============================================================

#' Tiled segmentation engine with seam-aware merging and global relabeling
#'
#' Runs a complete tiled segmentation workflow on a large `SpatRaster`:
#' (1) tile to disk with overlap, (2) segment each tile to disk with globally
#' unique IDs, (3) detect seam adjacencies, (4) compute per-segment means from
#' tiles, (5) merge adjacent seam segments by similarity, (6) apply merges
#' streamingly, and (7) relabel all segments globally to consecutive IDs.
#'
#' This function is intended as a robust, memory-safe orchestration layer for
#' classical OBIA workflows where full-scene segmentation is infeasible in RAM.
#'
#' @param x A `SpatRaster` to segment (single- or multi-band).
#' @param segment_fun A function implementing segmentation on a single tile.
#'   Must accept a `SpatRaster` as first argument and return a 1-layer integer
#'   segmentation raster.
#' @param seg_args A named list of additional arguments passed to `segment_fun`.
#' @param tile_size Integer. Tile size in pixels (rows/cols) excluding buffer.
#' @param buffer Integer. Overlap (in pixels) added around each tile.
#' @param seam_thr Numeric scalar. Threshold for merging seam-adjacent segments.
#'   Applied to Euclidean distance between per-segment mean feature vectors.
#' @param out_file Character. Output file path for the final segmentation.
#' @param tile_dir Character. Directory used for intermediate tile products.
#' @param cleanup_tiles Logical. If `TRUE`, remove the raw image tiles after
#'   segmentation.
#' @param cleanup_seg_tiles Logical. If `TRUE`, remove segmented tile files once
#'   the final output is produced.
#' @param verbose Logical. If `TRUE`, print progress messages.
#'
#' @return A file-backed `SpatRaster` pointing to `out_file`, containing a
#'   1-layer integer segmentation with globally consecutive IDs.
#'
#' @details
#'
#' See \code{vignette("tiled-segmentation-workflow")} for an end-to-end explanation
#' of the tiling and seam-merge framework.
#'
#' **Seam-aware merging.** Seam adjacencies are extracted only from the overlap
#' zone (buffer) of each tile to identify segment labels that touch across tile
#' borders. Segment means are computed from the original image values using
#' tile-wise streaming (no VRT reads). Adjacent segments are merged when their
#' mean vectors are within `seam_thr`, using transitive closure (union-find).
#'
#' **Fast path when no merges are needed.** If no seam adjacencies are found
#' (or if all candidate adjacencies are invalid due to missing means), the
#' function mosaics the segmented tiles directly to `out_file` without merging.
#'
#' **Streaming-safe output.** Before applying the merge map, the segmented tile
#' mosaic is materialized to a real file-backed raster to ensure stable block
#' reads/writes. Merging and final relabeling are performed using streaming
#' reads/writes to support large rasters.
#'
#' @seealso
#' [make_tiles_disk()], [segment_tiles()], `make_tile_seam_mask()`,
#' `extract_seam_pairs()`, [segment_means_from_tiles()],
#' [merge_by_adjacency()], [apply_merge_map()], [relabel_segments_global()]
#'
#' @export
#'
segmenter_tile_engine <- function(
    x,
    segment_fun,
    seg_args,
    tile_size = 2048,
    buffer = 64,
    seam_thr = 0.8,
    out_file = "final_seg.tif",
    tile_dir = tempdir(),
    cleanup_tiles = FALSE,
    cleanup_seg_tiles = FALSE,
    verbose=TRUE
) {
  stopifnot(inherits(x, "SpatRaster"))

  tile_dir <- normalizePath(tile_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(tile_dir)) dir.create(tile_dir, recursive = TRUE)


  tile_meta <- make_tiles_disk(x, tile_size = tile_size, buffer = buffer,
                               prefix = "tile", dir = tile_dir, overwrite = TRUE,
                               verbose=verbose)

  seg_dir <- file.path(tile_dir, "seg_tiles")

  if (!dir.exists(seg_dir)) dir.create(seg_dir, recursive = TRUE)

  seg_tiles <- segment_tiles(
    tile_meta,
    segment_fun,
    seg_args,
    out_dir = seg_dir,
    out_prefix = "seg",
    cleanup_tiles = cleanup_tiles,
    verbose=verbose
  )

  bad <- seg_tiles[!file.exists(seg_tiles)]
  if (length(bad)) stop("Missing segmented tile files: ", paste(bad, collapse = ", "))

  if(verbose) message("\n## [3/7] Making tile seam mask\n")

  # adjacency along seams
  adj_all <- NULL
  for (i in seq_along(seg_tiles)) {
    segt <- terra::rast(seg_tiles[i])
    seam_mask <- make_tile_seam_mask(segt, tile_meta[i, ], seam_width = buffer)
    adj <- extract_seam_pairs(segt, seam_mask = seam_mask, directions = 8)
    if (nrow(adj) > 0) adj_all <- if (is.null(adj_all)) adj else rbind(adj_all, adj)
  }



  # nothing to merge -> just mosaic
  if (is.null(adj_all) || !nrow(adj_all)) {

    if(verbose) message("\n## Nothing to merge! Mosaicking all files\n")

    terra::writeRaster(terra::vrt(seg_tiles), out_file, overwrite = TRUE,
                       wopt = list(datatype="INT4S"),
                       gdal = .default_gdal_opts())
    out <- terra::rast(out_file)
    if (cleanup_seg_tiles) try(file.remove(seg_tiles), silent = TRUE)
    return(out)
  }

  adj_all <- unique(adj_all)
  ids <- unique(c(adj_all))

  if(verbose) message("\n## [4/7] Computing segment means from tiles\n")

  # means computed from tiles (no VRT reads)
  means <- segment_means_from_tiles(x, seg_tiles, ids, block_nrows = 512)

  valid_ids <- rownames(means)[rowSums(is.na(means)) == 0]
  adj_all <- adj_all[adj_all[,1] %in% valid_ids & adj_all[,2] %in% valid_ids, , drop=FALSE]

  if (!nrow(adj_all)) {
    terra::writeRaster(terra::vrt(seg_tiles), out_file, overwrite = TRUE,
                       wopt = list(datatype="INT4S"),
                       gdal = .default_gdal_opts())
    out <- terra::rast(out_file)
    if (cleanup_seg_tiles) try(file.remove(seg_tiles), silent = TRUE)
    return(out)
  }

  if(verbose) message("\n## [5/7] Merging by adjacency\n")

  merge_map <- merge_by_adjacency(adj_all, means, seam_thr)

  # materialize mosaic to a real GeoTIFF before streaming merge map
  seg_all_file <- tempfile(fileext = ".tif")
  terra::writeRaster(terra::vrt(seg_tiles), seg_all_file, overwrite=TRUE,
                     wopt=list(datatype="INT4S"),
                     gdal=.default_gdal_opts())
  seg_all <- terra::rast(seg_all_file)

  if(verbose) message("\n## [6/7] Applying merge map\n")

  tmp_merged <- tempfile(fileext = ".tif")
  apply_merge_map(seg_all, merge_map, tmp_merged)

  if(verbose) message("\n## [7/7] Relabelling segments\n")

  out <- relabel_segments_global(terra::rast(tmp_merged), out_file)

  if (cleanup_seg_tiles) try(file.remove(seg_tiles), silent = TRUE)

  out
}

# ============================================================
# 8) User-facing wrappers
# ============================================================

#' Tiled Felzenszwalb-Huttenlocher (FH) segmentation with seam-aware merging
#'
#' Convenience wrapper around [segmenter_tile_engine()] that runs the package's
#' FH segmenter on disk-backed tiles and performs seam-aware merging and global
#' relabeling to produce a single, consistent segmentation raster.
#'
#' @param x A `SpatRaster` to segment.
#' @param tile_size Integer. Tile size in pixels (rows/cols) excluding buffer.
#' @param buffer Integer. Overlap (in pixels) added around each tile.
#' @param seam_thr Numeric. Threshold for merging seam-adjacent segments,
#'   applied to Euclidean distance between per-segment mean band vectors.
#' @param out_file Character. Output file path for the final segmentation.
#' @param tile_dir Character. Directory used for intermediate tile products.
#' @param cleanup_tiles Logical. If `TRUE`, remove raw image tile files after
#'   segmentation.
#' @param cleanup_seg_tiles Logical. If `TRUE`, remove segmented tile files once
#'   the final output is produced.
#'
#' @param k Numeric. FH scale parameter controlling the degree of merging
#'   (larger values generally yield larger segments).
#' @param min_size Integer. Minimum segment size (in pixels) enforced by FH.
#' @param eight Logical. If `TRUE`, use 8-neighborhood connectivity; otherwise
#'   use 4-neighborhood.
#' @param scale_bands Logical. If `TRUE`, standardize bands before segmentation.
#' @param smooth Integer. Optional spatial smoothing window size (in pixels);
#'   `0` disables smoothing.
#'
#' @return A file-backed `SpatRaster` pointing to `out_file`, containing a
#'   1-layer integer segmentation with globally consecutive IDs.
#'
#' @details
#' This function sets `segment_fun = fh_segmenter` and forwards FH parameters
#' through `seg_args` into [segmenter_tile_engine()].
#'
#' @seealso
#' [segmenter_tile_engine()], [fh_segmenter]
#'
#' @export
#'
fh_segmenter_tile <- function(x,
                              tile_size = 2048,
                              buffer = 64,
                              seam_thr = 0.7,
                              out_file = "fh_tiled_merged.tif",
                              tile_dir = tempdir(),
                              cleanup_tiles = FALSE,
                              cleanup_seg_tiles = FALSE,
                              k = 1.0,
                              min_size = 50,
                              eight = TRUE,
                              scale_bands = TRUE,
                              smooth = 0) {

  segmenter_tile_engine(
    x = x,
    segment_fun = fh_segmenter,
    seg_args = list(k=k, min_size=min_size, eight=eight,
                    scale_bands=scale_bands, smooth=smooth),
    tile_size = tile_size,
    buffer = buffer,
    seam_thr = seam_thr,
    out_file = out_file,
    tile_dir = tile_dir,
    cleanup_tiles = cleanup_tiles,
    cleanup_seg_tiles = cleanup_seg_tiles
  )
}

#' Tiled FH + region mean-shift segmentation with seam-aware merging
#'
#' Convenience wrapper around [segmenter_tile_engine()] that runs a two-stage
#' segmentation per tile: (1) Felzenszwalb-Huttenlocher (FH) superpixels and
#' (2) mean-shift clustering on region features, followed by seam-aware
#' merging and global relabeling.
#'
#' This wrapper provides tuned defaults for high-resolution imagery and is
#' intended as a robust, scalable OBIA-style segmenter for large rasters.
#'
#' @param x A `SpatRaster` to segment.
#' @param tile_size Integer. Tile size in pixels (rows/cols) excluding buffer.
#' @param buffer Integer. Overlap (in pixels) added around each tile.
#' @param seam_thr Numeric. Threshold for merging seam-adjacent segments,
#'   applied to Euclidean distance between per-segment mean band vectors.
#' @param out_file Character. Output file path for the final segmentation.
#' @param tile_dir Character. Directory used for intermediate tile products.
#' @param cleanup_tiles Logical. If `TRUE`, remove raw image tile files after
#'   segmentation.
#' @param cleanup_seg_tiles Logical. If `TRUE`, remove segmented tile files once
#'   the final output is produced.
#'
#' @param scale_bands Logical. If `TRUE`, standardize bands before processing.
#' @param smooth Integer. Optional smoothing window size (pixels) applied to
#'   the input prior to FH; `0` disables smoothing.
#'
#' @param fh_k Numeric. FH scale parameter (larger values generally yield
#'   larger base regions).
#' @param fh_min_size Integer. Minimum FH region size in pixels.
#' @param eight Logical. If `TRUE`, use 8-neighborhood connectivity; otherwise
#'   use 4-neighborhood.
#'
#' @param ms_dim Integer. Number of dimensions used for region feature
#'   vectors prior to mean-shift.
#' @param ms_ranger Numeric. Mean-shift range/bandwidth scaling parameter used
#'   by the region mean-shift stage.
#' @param ms_hs Numeric. Mean-shift spatial/kernel bandwidth (implementation-
#'   specific; forwarded to the underlying mean-shift routine).
#' @param ms_max_iter Integer. Maximum mean-shift iterations.
#' @param ms_eps Numeric. Convergence tolerance for mean-shift updates.
#' @param mode_merge Numeric. Post-processing threshold controlling merging of
#'   similar modes/clusters after mean-shift.
#' @param final_min_size Integer. Minimum final segment size (pixels) enforced
#'   after region clustering and mode merging.
#'
#' @return A file-backed `SpatRaster` pointing to `out_file`, containing a
#'   1-layer integer segmentation with globally consecutive IDs.
#'
#' @details
#' This function sets `segment_fun = fh_meanshift_segmenter` and forwards all
#' algorithm parameters through `seg_args` into [segmenter_tile_engine()].
#'
#' @seealso
#' [segmenter_tile_engine()], [fh_meanshift_segmenter]
#'
#' @export
#'
fh_meanshift_segmenter_tile <- function(x,
                                        tile_size = 1536,
                                        buffer = 96,
                                        seam_thr = 0.8,
                                        out_file = "fh_meanshift_tiled_merged.tif",
                                        tile_dir = tempdir(),
                                        cleanup_tiles = FALSE,
                                        cleanup_seg_tiles = FALSE,
                                        # --- preprocessing ---
                                        scale_bands = TRUE,
                                        smooth = 3,
                                        # --- FH ---
                                        fh_k = 0.5,
                                        fh_min_size = 20,
                                        eight = TRUE,
                                        # --- MeanShift on regions ---
                                        ms_dim = 3,
                                        ms_ranger = 0.15,
                                        ms_hs = 12,
                                        ms_max_iter = 10,
                                        ms_eps = 1e-3,
                                        mode_merge = 0.6,
                                        final_min_size = 80) {

  stopifnot(inherits(x, "SpatRaster"))

  segmenter_tile_engine(
    x = x,
    segment_fun = fh_meanshift_segmenter,
    seg_args = list(
      # preprocessing
      scale_bands = scale_bands,
      smooth = smooth,
      # FH
      fh_k = fh_k,
      fh_min_size = fh_min_size,
      eight = eight,
      # MeanShift on regions
      ms_dim = ms_dim,
      ms_ranger = ms_ranger,
      ms_hs = ms_hs,
      ms_max_iter = ms_max_iter,
      ms_eps = ms_eps,
      mode_merge = mode_merge,
      final_min_size = final_min_size
    ),
    tile_size = tile_size,
    buffer = buffer,
    seam_thr = seam_thr,
    out_file = out_file,
    tile_dir = tile_dir,
    cleanup_tiles = cleanup_tiles,
    cleanup_seg_tiles = cleanup_seg_tiles
  )
}
