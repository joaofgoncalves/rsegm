#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

namespace{

  static inline int qbin(double x, double lo, double hi, int nbins) {
    if (std::isnan(x)) return -1;
    if (hi <= lo) return 0;
    double t = (x - lo) / (hi - lo);
    if (t <= 0.0) return 0;
    if (t >= 1.0) return nbins - 1;
    int b = (int)std::floor(t * nbins);
    if (b < 0) b = 0;
    if (b >= nbins) b = nbins - 1;
    return b;
  }

  static inline void add_block_hist(std::vector<int>& out,
                                    const std::vector<int>& qimg,
                                    int nb, int nbins,
                                    int nrow, int ncol,
                                    int r0, int c0, int bs) {
    std::fill(out.begin(), out.end(), 0);
    for (int dr = 0; dr < bs; ++dr) {
      int r = r0 + dr;
      if (r < 0 || r >= nrow) continue;
      int base = r * ncol;
      for (int dc = 0; dc < bs; ++dc) {
        int c = c0 + dc;
        if (c < 0 || c >= ncol) continue;
        int idx = base + c;
        int off = idx * nb;
        for (int b = 0; b < nb; ++b) {
          int qb = qimg[off + b];
          if (qb >= 0) out[b * nbins + qb] += 1;
        }
      }
    }
  }

  // score = sum_i block[i] * ((hist_k[i]+alpha)/(size_k+alpha*H))
  static inline double block_score_prob(const std::vector<int>& block_hist,
                                        const std::vector<int>& hist_k,
                                        int size_k,
                                        double alpha) {
    const int H = (int)block_hist.size();
    const double denom = (double)size_k + alpha * (double)H;
    double s = 0.0;
    for (int i = 0; i < H; ++i) {
      int bc = block_hist[i];
      if (bc == 0) continue;
      double pk = ((double)hist_k[i] + alpha) / denom;
      s += (double)bc * pk;
    }
    return s;
  }

}

// [[Rcpp::export]]
IntegerVector seeds_segmenter_cpp(const NumericMatrix& img,
                                    int nrow, int ncol,
                                    int step,
                                    int nbins = 16,
                                    IntegerVector block_sizes = IntegerVector::create(8,4,2,1),
                                    int iters_per_level = 5,
                                    int boundary_samples = 64,
                                    double alpha = 1.0,
                                    bool verbose = false) {
  const int ncell = img.nrow();
  const int nb    = img.ncol();
  if (ncell != nrow * ncol) stop("img rows must equal nrow*ncol");
  if (step < 1) stop("step must be >= 1");
  if (nbins < 2) stop("nbins must be >= 2");

  // per-band min/max
  std::vector<double> lo(nb, R_PosInf), hi(nb, R_NegInf);
  for (int i = 0; i < ncell; ++i) {
    for (int b = 0; b < nb; ++b) {
      double x = img(i, b);
      if (std::isnan(x)) continue;
      if (x < lo[b]) lo[b] = x;
      if (x > hi[b]) hi[b] = x;
    }
  }
  for (int b = 0; b < nb; ++b) {
    if (!R_finite(lo[b]) || !R_finite(hi[b])) { lo[b] = 0.0; hi[b] = 1.0; }
    if (hi[b] <= lo[b]) hi[b] = lo[b] + 1.0;
  }

  // quantize
  std::vector<int> qimg((size_t)ncell * (size_t)nb);
  for (int i = 0; i < ncell; ++i) {
    for (int b = 0; b < nb; ++b) {
      qimg[(size_t)i * nb + b] = qbin(img(i, b), lo[b], hi[b], nbins);
    }
  }

  // initial grid labels
  const int gx = step;
  const int gy = step;
  const int nlabx = (ncol + gx - 1) / gx;
  const int nlaby = (nrow + gy - 1) / gy;
  const int K = nlabx * nlaby;

  std::vector<int> lab(ncell);
  for (int r = 0; r < nrow; ++r) {
    int ry = std::min(nlaby - 1, r / gy);
    for (int c = 0; c < ncol; ++c) {
      int cx = std::min(nlabx - 1, c / gx);
      lab[r * ncol + c] = ry * nlabx + cx; // 0..K-1
    }
  }

  const int H = nb * nbins;
  std::vector< std::vector<int> > hist(K, std::vector<int>(H, 0));
  std::vector<int> size(K, 0);

  for (int idx = 0; idx < ncell; ++idx) {
    int k = lab[idx];
    size[k] += 1;
    int off = idx * nb;
    for (int b = 0; b < nb; ++b) {
      int qb = qimg[off + b];
      if (qb >= 0) hist[k][b * nbins + qb] += 1;
    }
  }

  // deterministic RNG (CRAN-safe)
  uint32_t z = 2166136261u;
  auto rnd_idx = [&]() -> int {
    z = 1664525u * z + 1013904223u;
    return (int)(z % (uint32_t)ncell);
  };

  std::vector<int> block_hist(H);

  auto neighbor_label = [&](int idx, int kA) -> int {
    int r = idx / ncol;
    int c = idx - r * ncol;

    // check 4-neighbors, return first different label
    if (c + 1 < ncol) { int k = lab[idx + 1];     if (k != kA) return k; }
    if (c - 1 >= 0)   { int k = lab[idx - 1];     if (k != kA) return k; }
    if (r + 1 < nrow) { int k = lab[idx + ncol];  if (k != kA) return k; }
    if (r - 1 >= 0)   { int k = lab[idx - ncol];  if (k != kA) return k; }
    return -1;
  };

  for (int lvl = 0; lvl < block_sizes.size(); ++lvl) {
    const int bs = block_sizes[lvl];
    if (bs < 1) continue;

    for (int it = 0; it < iters_per_level; ++it) {
      if (verbose) Rcpp::checkUserInterrupt();

      // scale sampling with K (but cap for safety)
      const int nsamp = std::min(2000000, std::max(5000, boundary_samples * K));

      for (int s = 0; s < nsamp; ++s) {
        int idx = rnd_idx();
        int kA = lab[idx];
        int kB = neighbor_label(idx, kA);
        if (kB < 0) continue; // not a boundary sample

        int r = idx / ncol;
        int c = idx - r * ncol;

        // align block to bs grid around idx
        int r0 = (r / bs) * bs;
        int c0 = (c / bs) * bs;

        add_block_hist(block_hist, qimg, nb, nbins, nrow, ncol, r0, c0, bs);

        double sA = block_score_prob(block_hist, hist[kA], size[kA], alpha);
        double sB = block_score_prob(block_hist, hist[kB], size[kB], alpha);

        if (sB <= sA) continue;

        // accept move: only pixels currently in kA within block change to kB
        for (int dr = 0; dr < bs; ++dr) {
          int rr = r0 + dr;
          if (rr < 0 || rr >= nrow) continue;
          int base = rr * ncol;
          for (int dc = 0; dc < bs; ++dc) {
            int cc = c0 + dc;
            if (cc < 0 || cc >= ncol) continue;
            int p = base + cc;
            if (lab[p] != kA) continue;

            lab[p] = kB;
            size[kA] -= 1;
            size[kB] += 1;

            int off = p * nb;
            for (int b = 0; b < nb; ++b) {
              int qb = qimg[off + b];
              if (qb >= 0) {
                hist[kA][b * nbins + qb] -= 1;
                hist[kB][b * nbins + qb] += 1;
              }
            }
          }
        }
      }
    }
  }

  IntegerVector out(ncell);
  for (int i = 0; i < ncell; ++i) out[i] = lab[i] + 1;
  return out;
}
