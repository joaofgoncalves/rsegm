> ⚠️ **Warning: development version**
>
> This R package is currently in active development.\
> Features may be incomplete, change without notice, or behave unexpectedly.\
> Use with caution, and do not rely on it for production workflows.

# Geospatial image segmentation in R

This package provides scalable, high-performance tools for **geospatial image segmentation** built on top of **terra**, **Rcpp**, and modern block-wise / tiled processing strategies. It is designed for object-based image analysis (OBIA) workflows on large raster datasets, including satellite and aerial imagery.

The core focus is on producing **spatially coherent, integer-labeled segment rasters** that integrate cleanly with downstream raster and vector analysis in R.

![](man/figures/showcase.png)

## Main features

-   **Graph-based and hybrid segmentation algorithms**
    -   Fast Felzenszwalb–Huttenlocher (FH) graph-based segmentation
    -   Hybrid FH + region-level Mean Shift segmentation for higher-quality objects
-   **Large-raster support**
    -   Robust tiled processing for rasters that do not fit comfortably in memory
    -   Seam-aware merging to avoid boundary artifacts
-   **Efficient C++ backends**
    -   Computationally intensive steps implemented in Rcpp for speed and scalability
-   **terra-native design**
    -   Uses `SpatRaster` throughout
    -   Produces standard GeoTIFF outputs with integer segment IDs

## Segmentation methods

The package currently includes:

-   **FH graph-based segmentation**\
    A fast, noise-robust method often used to generate superpixels or initial regions, suitable for large multi-band rasters.

-   **FH + mean shift hybrid segmentation**\
    A two-stage approach combining FH over-segmentation with region-level Mean Shift refinement to merge spectrally similar regions efficiently, producing object-like segments well suited for OBIA workflows.

All methods return a single-layer raster where each cell contains an integer segment identifier.

**Example**

``` r
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

For large scenes, segmentation can be executed in a **tiled workflow** that:

1.  Splits the input raster into buffered tiles on disk\
2.  Segments each tile independently\
3.  Ensures globally unique segment IDs across tiles\
4.  Detects and merges adjacent segments across tile seams\
5.  Produces a single, seamless output raster

This approach allows segmentation of very large rasters while maintaining correctness and reproducibility.

## Typical workflow with tiling

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

## Tiling and segmentation parameters

This section explains the main parameters controlling tiled segmentation and provides practical guidance for selecting appropriate values.

### `tile_size`

Defines the **spatial size (in pixels)** of each processing tile.\
Each tile is segmented independently and later merged.

-   Larger values reduce the number of seams but require more memory
-   Smaller values lower memory usage but increase the number of tile boundaries

This parameter controls the **fundamental processing unit** of the tiled workflow.

------------------------------------------------------------------------

### `buffer`

Defines the **overlap (in pixels)** added around each tile.

-   Prevents artificial boundaries at tile edges
-   Ensures sufficient spatial context for segments near tile borders
-   Used during both segmentation and seam reconciliation

The buffer should generally be **larger than the expected object size** and larger than the segmentation neighborhood.

------------------------------------------------------------------------

### `k`

Segmentation **scale parameter** for graph-based segmentation.

-   Lower values produce many small segments (over-segmentation)
-   Higher values produce fewer, larger segments

While `k` is independent of tiling, it strongly influences how visible tile seams become if set too low.

### Recommended settings

| Scenario | tile_size | buffer | k | Notes |
|----|----|----|----|----|
| Small raster (\< 5k x 5k) | 2048–4096 | 64 | 0.4–0.7 | Larger tiles reduce seam handling |
| Large raster, limited RAM | 1024–2048 | 64–96 | 0.6–1.0 | Balance memory use and seam robustness |
| High-resolution imagery (\<=0.5 m) | 2048 | 96–128 | 0.3–0.6 | Larger buffer needed for fine detail |
| Coarse-resolution imagery (\>= 10 m) | 4096 | 32–64 | 0.8–1.5 | Fewer seams, larger objects |
| Urban / textured scenes | 2048 | \>= 96 | 0.4–0.7 | Buffer critical to avoid edge artifacts |
| Homogeneous landscapes | 4096 | 32–64 | 0.8–1.2 | Large tiles preferred |

### Practical guidelines

-   Increase **`tile_size` first** if tiling patterns are visible
-   Increase **`buffer`** if segment boundaries align with tile edges
-   Adjust **`k` last** to control object size once seams are handled
-   Ensure:\
    **buffer \>= expected object radius**\
    **tile_size \>\> buffer**

The result is a GeoTIFF with globally consistent segment IDs, suitable for visualization, statistics, or conversion to vector objects.

## Design principles

-   Favor **explicitness over magic**: segment IDs, tiling, and merging steps are transparent.
-   Be **robust to raster size and storage mode** (in-memory vs file-backed).
-   Keep outputs **simple and interoperable**: integer labels, one layer, no hidden state.

## Dependencies

-   **terra** for raster I/O and spatial operations
-   **Rcpp** for performance-critical algorithms
-   Base R (stats, utils)

## License

See the `DESCRIPTION` file for licensing and package metadata.
