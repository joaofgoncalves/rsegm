// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <limits>
#include <tuple>
#include <utility>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

namespace{

  // Segment structure to hold region information
  struct Segment {
    int id;
    int n_pixels;
    arma::vec mean_values;      // Mean spectral values for each band
    arma::vec sum_values;       // Sum of spectral values
    arma::vec sum_sq_values;    // Sum of squared values for variance
    std::unordered_set<int> pixels;  // Pixel indices in this segment
    std::unordered_set<int> neighbors; // Neighboring segment IDs
    double perimeter;
    double bbox_perimeter;

    Segment(int id_, int nbands) :
      id(id_), n_pixels(0), perimeter(0.0), bbox_perimeter(0.0) {
      mean_values.zeros(nbands);
      sum_values.zeros(nbands);
      sum_sq_values.zeros(nbands);
    }
  };

  // Merge candidate structure for priority queue
  struct MergeCandidate {
    int seg1_id;
    int seg2_id;
    double fusion_factor;

    // Greater than for min-heap (lowest fusion factor has highest priority)
    bool operator>(const MergeCandidate& other) const {
      return fusion_factor > other.fusion_factor;
    }
  };

  // Calculate perimeter of a segment
  double calculatePerimeter(const std::unordered_set<int>& pixels,
                            int nrows, int ncols) {
    double perim = 0.0;
    for (int idx : pixels) {
      int r = idx / ncols;
      int c = idx % ncols;

      // Check 4-connectivity neighbors
      std::vector<std::pair<int,int>> neighbors = {
        {r-1, c}, {r+1, c}, {r, c-1}, {r, c+1}
      };

      for (const auto& nb : neighbors) {
        int nr = nb.first;
        int nc = nb.second;
        if (nr < 0 || nr >= nrows || nc < 0 || nc >= ncols) {
          perim += 1.0;  // Edge of image
        } else {
          int nb_idx = nr * ncols + nc;
          if (pixels.find(nb_idx) == pixels.end()) {
            perim += 1.0;  // Neighbor not in segment
          }
        }
      }
    }
    return perim;
  }

  // Calculate bounding box perimeter
  double calculateBBoxPerimeter(const std::unordered_set<int>& pixels, int ncols) {
    if (pixels.empty()) return 0.0;

    int min_r = std::numeric_limits<int>::max();
    int max_r = std::numeric_limits<int>::min();
    int min_c = std::numeric_limits<int>::max();
    int max_c = std::numeric_limits<int>::min();

    for (int idx : pixels) {
      int r = idx / ncols;
      int c = idx % ncols;
      min_r = std::min(min_r, r);
      max_r = std::max(max_r, r);
      min_c = std::min(min_c, c);
      max_c = std::max(max_c, c);
    }

    int width = max_c - min_c + 1;
    int height = max_r - min_r + 1;
    return 2.0 * (width + height);
  }

  // Calculate spectral heterogeneity increase from merging
  double calcSpectralHeterogeneity(const Segment& s1, const Segment& s2,
                                   const arma::vec& band_weights) {
    int n1 = s1.n_pixels;
    int n2 = s2.n_pixels;
    int n3 = n1 + n2;

    if (n3 == 0) return 0.0;

    double h_color = 0.0;
    int nbands = s1.mean_values.n_elem;

    for (int b = 0; b < nbands; ++b) {
      double var1 = (n1 > 1) ?
      (s1.sum_sq_values(b) / n1 - std::pow(s1.mean_values(b), 2)) : 0.0;
      double var2 = (n2 > 1) ?
      (s2.sum_sq_values(b) / n2 - std::pow(s2.mean_values(b), 2)) : 0.0;

      // Merged variance
      double sum3 = s1.sum_values(b) + s2.sum_values(b);
      double sum_sq3 = s1.sum_sq_values(b) + s2.sum_sq_values(b);
      double mean3 = sum3 / n3;
      double var3 = (sum_sq3 / n3 - std::pow(mean3, 2));

      // Standard deviation weighted by number of pixels
      double std1 = std::sqrt(std::max(0.0, var1));
      double std2 = std::sqrt(std::max(0.0, var2));
      double std3 = std::sqrt(std::max(0.0, var3));

      h_color += band_weights(b) * (n3 * std3 - n1 * std1 - n2 * std2);
    }

    return h_color;
  }

  // Calculate shape heterogeneity increase from merging
  double calcShapeHeterogeneity(const Segment& s1, const Segment& s2,
                                double l1, double l2, double l3,
                                double b1, double b2, double b3,
                                double w_compactness) {
    int n1 = s1.n_pixels;
    int n2 = s2.n_pixels;
    int n3 = n1 + n2;

    if (n3 == 0) return 0.0;

    // Compactness component
    double cmpct1 = (n1 > 0) ? l1 / std::sqrt(static_cast<double>(n1)) : 0.0;
    double cmpct2 = (n2 > 0) ? l2 / std::sqrt(static_cast<double>(n2)) : 0.0;
    double cmpct3 = (n3 > 0) ? l3 / std::sqrt(static_cast<double>(n3)) : 0.0;
    double h_cmpct = n3 * cmpct3 - n1 * cmpct1 - n2 * cmpct2;

    // Smoothness component
    double smooth1 = (n1 > 0 && b1 > 0) ? l1 / b1 : 0.0;
    double smooth2 = (n2 > 0 && b2 > 0) ? l2 / b2 : 0.0;
    double smooth3 = (n3 > 0 && b3 > 0) ? l3 / b3 : 0.0;
    double h_smooth = n3 * smooth3 - n1 * smooth1 - n2 * smooth2;

    return w_compactness * h_cmpct + (1.0 - w_compactness) * h_smooth;
  }

  // Calculate fusion factor between two segments
  double calcFusionFactor(const Segment& s1, const Segment& s2,
                          const arma::vec& band_weights,
                          double w_color, double w_compactness,
                          int nrows, int ncols) {
    // Calculate merged properties
    std::unordered_set<int> merged_pixels = s1.pixels;
    merged_pixels.insert(s2.pixels.begin(), s2.pixels.end());

    double l1 = s1.perimeter;
    double l2 = s2.perimeter;
    double l3 = calculatePerimeter(merged_pixels, nrows, ncols);

    double b1 = s1.bbox_perimeter;
    double b2 = s2.bbox_perimeter;
    double b3 = calculateBBoxPerimeter(merged_pixels, ncols);

    // Spectral heterogeneity
    double h_color = calcSpectralHeterogeneity(s1, s2, band_weights);

    // Shape heterogeneity
    double h_shape = calcShapeHeterogeneity(s1, s2, l1, l2, l3,
                                            b1, b2, b3, w_compactness);

    // Fusion factor
    return w_color * h_color + (1.0 - w_color) * h_shape;
  }

  // Merge two segments
  void mergeSegments(Segment& target, const Segment& source,
                     int nrows, int ncols) {
    // Update pixel count
    target.n_pixels += source.n_pixels;

    // Update spectral statistics
    target.sum_values += source.sum_values;
    target.sum_sq_values += source.sum_sq_values;
    target.mean_values = target.sum_values / target.n_pixels;

    // Merge pixels
    target.pixels.insert(source.pixels.begin(), source.pixels.end());

    // Merge neighbors (excluding the merged segment)
    for (int nb : source.neighbors) {
      if (nb != target.id) {
        target.neighbors.insert(nb);
      }
    }
    target.neighbors.erase(source.id);

    // Recalculate geometric properties
    target.perimeter = calculatePerimeter(target.pixels, nrows, ncols);
    target.bbox_perimeter = calculateBBoxPerimeter(target.pixels, ncols);
  }

}

// Baatz-Schape Multi-Resolution Segmentation
// [[Rcpp::export]]
Rcpp::IntegerVector baatz_segmenter_cpp(
     const arma::mat& image,
     int nrows,
     int ncols,
     double scale_param,
     double color_weight = 0.9,
     double compactness_weight = 0.5,
     Rcpp::Nullable<Rcpp::NumericVector> band_weights = R_NilValue
 ) {
   int npixels = nrows * ncols;
   int nbands = image.n_cols;

   // Setup band weights
   arma::vec bw(nbands);
   if (band_weights.isNotNull()) {
     NumericVector bw_r(band_weights);
     for (int i = 0; i < nbands; ++i) {
       bw(i) = bw_r(i);
     }
   } else {
     bw.fill(1.0);
   }

   Rcpp::Rcout << "Initializing segmentation..." << std::endl;
   Rcpp::Rcout << "  Image: " << nrows << "x" << ncols << " pixels, "
               << nbands << " bands" << std::endl;
   Rcpp::Rcout << "  Scale parameter: " << scale_param << std::endl;
   Rcpp::Rcout << "  Color weight: " << color_weight << std::endl;

   // Initialize: each pixel is a segment
   std::unordered_map<int, Segment> segments;
   segments.reserve(npixels);  // Reserve space for efficiency
   IntegerVector labels(npixels);

   for (int i = 0; i < npixels; ++i) {
     // Use emplace to construct Segment directly in the map
     // This avoids the need for a default constructor
     auto result = segments.emplace(std::piecewise_construct,
                                    std::forward_as_tuple(i),
                                    std::forward_as_tuple(i, nbands));

     Segment& seg = result.first->second;
     seg.n_pixels = 1;
     seg.pixels.insert(i);

     // Initialize spectral values
     for (int b = 0; b < nbands; ++b) {
       double val = image(i, b);
       seg.sum_values(b) = val;
       seg.sum_sq_values(b) = val * val;
       seg.mean_values(b) = val;
     }

     // Calculate initial perimeter and bbox
     seg.perimeter = 4.0;  // Single pixel has perimeter 4
     seg.bbox_perimeter = 4.0;

     labels(i) = i;
   }

   // Build initial neighbor relationships
   Rcpp::Rcout << "Building neighbor relationships..." << std::endl;
   for (int r = 0; r < nrows; ++r) {
     for (int c = 0; c < ncols; ++c) {
       int idx = r * ncols + c;

       // 4-connectivity
       std::vector<std::pair<int,int>> nbs = {
         {r-1, c}, {r+1, c}, {r, c-1}, {r, c+1}
       };

       for (const auto& nb : nbs) {
         if (nb.first >= 0 && nb.first < nrows &&
             nb.second >= 0 && nb.second < ncols) {
           int nb_idx = nb.first * ncols + nb.second;
           segments.at(idx).neighbors.insert(nb_idx);
         }
       }
     }
   }

   // Iterative merging
   double scale_threshold = scale_param * scale_param;
   int iteration = 0;
   bool merged = true;

   Rcpp::Rcout << "Starting region merging iterations..." << std::endl;

   while (merged) {
     merged = false;
     iteration++;

     std::priority_queue<MergeCandidate,
                         std::vector<MergeCandidate>,
                         std::greater<MergeCandidate>> merge_queue;

     std::unordered_set<int> processed;

     // Find best merge candidate for each segment
     for (auto& seg_pair : segments) {
       Segment& seg = seg_pair.second;

       if (processed.find(seg.id) != processed.end()) continue;

       double min_ff = std::numeric_limits<double>::max();
       int best_neighbor = -1;

       for (int nb_id : seg.neighbors) {
         if (segments.find(nb_id) == segments.end()) continue;
         if (processed.find(nb_id) != processed.end()) continue;

         const Segment& nb = segments.at(nb_id);
         double ff = calcFusionFactor(seg, nb, bw, color_weight,
                                      compactness_weight, nrows, ncols);

         if (ff < min_ff) {
           min_ff = ff;
           best_neighbor = nb_id;
         }
       }

       if (best_neighbor >= 0 && min_ff < scale_threshold) {
         merge_queue.push({seg.id, best_neighbor, min_ff});
       }
     }

     // Process merges
     int merge_count = 0;
     while (!merge_queue.empty()) {
       MergeCandidate mc = merge_queue.top();
       merge_queue.pop();

       // Check if both segments still exist
       if (segments.find(mc.seg1_id) == segments.end() ||
           segments.find(mc.seg2_id) == segments.end()) {
         continue;
       }

       // Perform merge
       Segment& seg1 = segments.at(mc.seg1_id);
       Segment& seg2 = segments.at(mc.seg2_id);

       // Merge seg2 into seg1
       mergeSegments(seg1, seg2, nrows, ncols);

       // Update neighbor relationships
       for (int nb_id : seg2.neighbors) {
         auto it = segments.find(nb_id);
         if (it != segments.end() && nb_id != mc.seg1_id) {
           it->second.neighbors.erase(mc.seg2_id);
           it->second.neighbors.insert(mc.seg1_id);
         }
       }

       // Update labels
       for (int pix_id : seg2.pixels) {
         labels(pix_id) = mc.seg1_id;
       }

       // Remove seg2
       segments.erase(mc.seg2_id);

       merge_count++;
       merged = true;
     }

     if (iteration % 5 == 0 || !merged) {
       Rcpp::Rcout << "  Iteration " << iteration
                   << ": " << segments.size() << " segments, "
                   << merge_count << " merges" << std::endl;
     }

     if (segments.size() <= 1) break;
   }

   Rcpp::Rcout << "Segmentation complete!" << std::endl;
   Rcpp::Rcout << "  Final: " << segments.size() << " segments" << std::endl;

   // Relabel segments sequentially
   std::unordered_map<int, int> label_map;
   int new_label = 1;
   for (const auto& seg_pair : segments) {
     label_map[seg_pair.first] = new_label++;
   }

   for (int i = 0; i < npixels; ++i) {
     labels(i) = label_map[labels(i)];
   }

   return labels;
 }
