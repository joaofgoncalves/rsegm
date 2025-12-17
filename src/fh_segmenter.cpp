#include <Rcpp.h>
#include <R.h>     // ISNA, R_finite
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Rcpp;

namespace{

  struct Edge {
    int a;
    int b;
    float w;
  };

  struct DSU {
    std::vector<int> parent;
    std::vector<int> rank;
    std::vector<int> size;
    std::vector<float> intr; // "internal difference" per component root

    explicit DSU(int n)
      : parent(n), rank(n, 0), size(n, 1), intr(n, 0.0f) {
      for (int i = 0; i < n; ++i) parent[i] = i;
    }

    inline int find(int x) {
      // Robust, no path-halving (avoids parent[parent[x]] access)
      const int n = (int)parent.size();
      if (x < 0 || x >= n) Rcpp::stop("DSU::find: x out of bounds");

      int r = x;
      while (true) {
        int pr = parent[r];
        if (pr < 0 || pr >= n) Rcpp::stop("DSU::find: parent out of bounds");
        if (pr == r) break;
        r = pr;
      }

      // compress
      while (parent[x] != x) {
        int px = parent[x];
        parent[x] = r;
        x = px;
        if (x < 0 || x >= n) Rcpp::stop("DSU::find: compress step OOB");
      }

      return r;
    }

    inline void unite_roots(int ra, int rb, float w) {
      // ra and rb MUST be roots
      if (ra == rb) return;

      const int n = (int)parent.size();
      if (ra < 0 || ra >= n || rb < 0 || rb >= n) Rcpp::stop("DSU::unite: root OOB");
      if (parent[ra] != ra || parent[rb] != rb) Rcpp::stop("DSU::unite: non-root passed");

      if (rank[ra] < rank[rb]) std::swap(ra, rb);

      parent[rb] = ra;
      size[ra] += size[rb];

      // FH: internal difference of merged component = max(intr(A), intr(B), w)
      float m = intr[ra];
      if (intr[rb] > m) m = intr[rb];
      if (w > m) m = w;
      intr[ra] = m;

      if (rank[ra] == rank[rb]) rank[ra] += 1;
    }
  };

  static inline float tau(float internal, int comp_size, float k) {
    // FH threshold function
    return internal + k / (float)comp_size;
  }

  static inline bool pixel_valid_allbands(const double* img, int p, int nb, int npix) {
    for (int b = 0; b < nb; ++b) {
      const double v = img[p + npix * b]; // band-stacked (column-major by band)
      if (ISNA(v) || !R_finite(v)) return false;
    }
    return true;
  }

  static inline float pixdist(const double* img, int p, int q, int nb, int npix) {
    double acc = 0.0;
    for (int b = 0; b < nb; ++b) {
      const double d = img[p + npix * b] - img[q + npix * b];
      acc += d * d;
    }
    return (float)std::sqrt(acc);
  }
}


// [[Rcpp::export]]
IntegerVector fh_segmenter_cpp(NumericVector img,
                             int nrow, int ncol, int nb,
                             double k = 0.8,
                             int min_size = 30,
                             bool eight = true) {

  if (nrow <= 0 || ncol <= 0 || nb <= 0) Rcpp::stop("nrow/ncol/nb must be > 0");
  if (!(k > 0.0)) Rcpp::stop("k must be > 0");
  if (min_size < 1) min_size = 1;

  const size_t npix_sz = (size_t)nrow * (size_t)ncol;
  if (npix_sz > (size_t)INT_MAX) Rcpp::stop("image too large (npix > INT_MAX)");
  const int npix = (int)npix_sz;

  const size_t need = (size_t)npix * (size_t)nb;
  if ((size_t)img.size() != need) Rcpp::stop("img length must be nrow*ncol*nb");

  const double* pimg = REAL(img);
  const float kf = (float)k;

  // Valid mask: pixel is valid only if all bands are finite and non-NA
  std::vector<unsigned char> valid(npix, 0);
  for (int p = 0; p < npix; ++p) {
    valid[p] = pixel_valid_allbands(pimg, p, nb, npix) ? 1u : 0u;
  }

  // Build adjacency edges (only between valid pixels)
  std::vector<Edge> edges;
  // Reserve a safe upper bound (avoid pathological reserve)
  const size_t max_edges = (size_t)npix * (eight ? 4u : 2u);
  edges.reserve(max_edges);

  auto id = [ncol](int r, int c) { return r * ncol + c; };

  for (int r = 0; r < nrow; ++r) {
    for (int c = 0; c < ncol; ++c) {
      const int a = id(r, c);
      if (!valid[a]) continue;

      // right
      if (c + 1 < ncol) {
        const int b = id(r, c + 1);
        if (valid[b]) edges.push_back({a, b, pixdist(pimg, a, b, nb, npix)});
      }
      // down
      if (r + 1 < nrow) {
        const int b = id(r + 1, c);
        if (valid[b]) edges.push_back({a, b, pixdist(pimg, a, b, nb, npix)});
      }

      if (eight) {
        // down-right
        if (r + 1 < nrow && c + 1 < ncol) {
          const int b = id(r + 1, c + 1);
          if (valid[b]) edges.push_back({a, b, pixdist(pimg, a, b, nb, npix)});
        }
        // down-left
        if (r + 1 < nrow && c - 1 >= 0) {
          const int b = id(r + 1, c - 1);
          if (valid[b]) edges.push_back({a, b, pixdist(pimg, a, b, nb, npix)});
        }
      }
    }
  }

  // Sort edges by weight
  std::sort(edges.begin(), edges.end(),
            [](const Edge& e1, const Edge& e2) { return e1.w < e2.w; });

  DSU dsu(npix);

  // Mark invalid pixels as size=0 components so they never merge
  for (int p = 0; p < npix; ++p) {
    if (!valid[p]) dsu.size[p] = 0;
  }

  // PASS 1: FH criterion
  for (size_t i = 0; i < edges.size(); ++i) {
    const Edge& e = edges[i];

    // periodic interrupt check (keeps R responsive for large images)
    if ((i & 0xFFFF) == 0) Rcpp::checkUserInterrupt();

    int ra = dsu.find(e.a);
    int rb = dsu.find(e.b);
    if (ra == rb) continue;

    if (dsu.size[ra] == 0 || dsu.size[rb] == 0) continue;

    const float ta = tau(dsu.intr[ra], dsu.size[ra], kf);
    const float tb = tau(dsu.intr[rb], dsu.size[rb], kf);

    if (e.w <= (ta < tb ? ta : tb)) {
      dsu.unite_roots(ra, rb, e.w);
    }
  }

  // PASS 2: enforce minimum component size (standard FH post-process)
  if (min_size > 1) {
    for (size_t i = 0; i < edges.size(); ++i) {
      const Edge& e = edges[i];
      if ((i & 0xFFFF) == 0) Rcpp::checkUserInterrupt();

      int ra = dsu.find(e.a);
      int rb = dsu.find(e.b);
      if (ra == rb) continue;

      if (dsu.size[ra] == 0 || dsu.size[rb] == 0) continue;

      if (dsu.size[ra] < min_size || dsu.size[rb] < min_size) {
        dsu.unite_roots(ra, rb, e.w);
      }
    }
  }

  // Relabel (vector remap, stable, no hashing)
  IntegerVector out(npix, 0);
  std::vector<int> remap(npix, 0);
  int nextLab = 1;

  for (int p = 0; p < npix; ++p) {
    if (!valid[p]) { out[p] = 0; continue; }
    const int rp = dsu.find(p);
    int &lab = remap[rp];
    if (lab == 0) lab = nextLab++;
    out[p] = lab;
  }

  return out;
}
