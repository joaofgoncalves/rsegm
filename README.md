---
editor_options: 
  markdown: 
    wrap: 72
---

# Geospatial image segmentation in R

This package provides scalable, high-performance tools for **geospatial
image segmentation** built on top of **terra**, **Rcpp**, and modern
block-wise / tiled processing strategies. It is designed for
object-based image analysis (OBIA) workflows on large raster datasets,
including satellite and aerial imagery.

The core focus is on producing **spatially coherent, integer-labeled
segment rasters** that integrate cleanly with downstream raster and
vector analysis in R.

## Main features

-   **Graph-based and hybrid segmentation algorithms**
    -   Fast Felzenszwalb-Huttenlocher (FH) graph-based segmentation
    -   Hybrid FH + region-level Mean Shift segmentation for
        higher-quality objects
-   **Large-raster support**
    -   Robust tiled processing for rasters that do not fit comfortably
        in memory
    -   Seam-aware merging to avoid boundary artifacts
-   **Efficient C++ backends**
    -   Computationally intensive steps implemented in Rcpp for speed
        and scalability
-   **terra-native design**
    -   Uses `SpatRaster` throughout
    -   Produces standard GeoTIFF outputs with integer segment IDs

## Segmentation methods

The package currently includes:

-   **FH graph-based segmentation**\
    A fast, noise-robust method often used to generate superpixels or
    initial regions, suitable for large multi-band rasters.

-   **FH + mean shift hybrid segmentation**\
    A two-stage approach combining FH over-segmentation with
    region-level Mean Shift refinement to merge spectrally similar
    regions efficiently, producing object-like segments well suited for
    OBIA workflows.

All methods return a single-layer raster where each cell contains an
integer segment identifier.

**Example**

``` R
library(terra)
library(rsegm)

r <- rast("PlanetSAT_10m_Finland_Helsinki_UTM35.tif") 

print(r)

# Segment image using FH method
#
seg <- fh_segmenter(
  r,
  k = 0.5,
  min_size = 40,
  eight = TRUE,
  scale_bands = TRUE,
  smooth = 3
)

seg_pol   <- as.polygons(seg, dissolve = TRUE)
seg_lines <- as.lines(seg_pol)

rgb <- r[[1:3]]
rgb <- stretch(rgb)

plotRGB(rgb, r=1, g=2, b=3, stretch="lin")
plot(seg_lines, add=TRUE, col="yellow", lwd=0.8)
```

## Tiled processing [TBD]

For large scenes, segmentation can be executed in a **tiled workflow**
that:

1.  Splits the input raster into buffered tiles on disk\
2.  Segments each tile independently\
3.  Ensures globally unique segment IDs across tiles\
4.  Detects and merges adjacent segments across tile seams\
5.  Produces a single, seamless output raster

This approach allows segmentation of very large rasters while
maintaining correctness and reproducibility.

## Typical workflow

``` r

library(terra)
library(rsegm)

r <- rast("multispectral_image.tif")

seg <- fh_segmenter_tile(
  r,
  tile_size = 2048,
  buffer = 64,
  k = 0.7,
  min_size = 60,
  smooth = 3,
  out_file = "segmented.tif"
)
```

The result is a GeoTIFF with globally consistent segment IDs, suitable
for visualization, statistics, or conversion to vector objects.

## Design principles

-   Favor **explicitness over magic**: segment IDs, tiling, and merging
    steps are transparent.
-   Be **robust to raster size and storage mode** (in-memory vs
    file-backed).
-   Keep outputs **simple and interoperable**: integer labels, one
    layer, no hidden state.

## Dependencies

-   **terra** for raster I/O and spatial operations
-   **Rcpp** for performance-critical algorithms
-   Base R (stats, utils)

## License

See the `DESCRIPTION` file for licensing and package metadata.
