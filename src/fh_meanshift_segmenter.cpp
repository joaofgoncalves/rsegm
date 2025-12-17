// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <cstdint>

using namespace Rcpp;

namespace{

  // -----------------------------
  // DSU (generic)
  // -----------------------------
  struct DSU {
    std::vector<int> p, r, sz;
    DSU(int n): p(n), r(n,0), sz(n,1) { for(int i=0;i<n;++i) p[i]=i; }
    inline int find(int x){
      while(p[x]!=x){ p[x]=p[p[x]]; x=p[x]; }
      return x;
    }
    inline int unite(int a,int b){
      a=find(a); b=find(b);
      if(a==b) return a;
      if(r[a]<r[b]) std::swap(a,b);
      p[b]=a; sz[a]+=sz[b];
      if(r[a]==r[b]) r[a]++;
      return a;
    }
  };

  // -----------------------------
  // FH (Felzenszwalb-Huttenlocher) DSU
  // -----------------------------
  struct Edge { int a, b; float w; };

  struct DSUFH {
    std::vector<int> p, r, sz;
    std::vector<float> internal;
    DSUFH(int n): p(n), r(n,0), sz(n,1), internal(n,0.0f) { for(int i=0;i<n;++i) p[i]=i; }
    inline int find(int x){
      while(p[x]!=x){ p[x]=p[p[x]]; x=p[x]; }
      return x;
    }
    inline void unite(int a,int b,float w){
      a=find(a); b=find(b);
      if(a==b) return;
      if(r[a]<r[b]) std::swap(a,b);
      p[b]=a;
      sz[a]+=sz[b];
      internal[a] = std::max(std::max(internal[a], internal[b]), w);
      if(r[a]==r[b]) r[a]++;
    }
  };

  inline float tau(float internal, int size, float k){
    return internal + k / (float)size;
  }

  inline float pixDist(const float* img, int npix, int nb, int p, int q){
    double acc=0.0;
    for(int k=0;k<nb;++k){
      double d = (double)img[p + npix*k] - (double)img[q + npix*k];
      acc += d*d;
    }
    return (float)std::sqrt(acc);
  }

  inline std::uint64_t packEdge(int a, int b){
    if(a>b) std::swap(a,b);
    return ( (std::uint64_t)(std::uint32_t)a << 32 ) | (std::uint64_t)(std::uint32_t)b;
  }

  // -----------------------------
  // Key for grid hashing (dim <= 5)
  // -----------------------------
  struct Key {
    int x0,x1,x2,x3,x4,dim;
    bool operator==(Key const& o) const {
      if(dim!=o.dim) return false;
      if(x0!=o.x0 || x1!=o.x1) return false;
      if(dim<=2) return true;
      if(x2!=o.x2) return false;
      if(dim<=3) return true;
      if(x3!=o.x3) return false;
      if(dim<=4) return true;
      return x4==o.x4;
    }
  };

  struct KeyHash {
    std::size_t operator()(Key const& k) const noexcept {
      std::size_t h = 1469598103934665603ULL;
      auto mix = [&](int v){
        h ^= (std::size_t)(std::uint32_t)v + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
      };
      mix(k.dim); mix(k.x0); mix(k.x1);
      if(k.dim>2) mix(k.x2);
      if(k.dim>3) mix(k.x3);
      if(k.dim>4) mix(k.x4);
      return h;
    }
  };

  inline Key makeKey(const std::vector<float>& z, int i, int dim){
    Key k; k.dim=dim;
    k.x0 = (int)std::floor(z[(size_t)i*dim + 0]);
    k.x1 = (int)std::floor(z[(size_t)i*dim + 1]);
    k.x2 = (dim>2)? (int)std::floor(z[(size_t)i*dim + 2]) : 0;
    k.x3 = (dim>3)? (int)std::floor(z[(size_t)i*dim + 3]) : 0;
    k.x4 = (dim>4)? (int)std::floor(z[(size_t)i*dim + 4]) : 0;
    return k;
  }

  static void offsets_rec(int dim, int d, std::vector<int>& cur, std::vector< std::vector<int> >& out){
    if(d==dim){ out.push_back(cur); return; }
    for(int v=-1; v<=1; ++v){
      cur[d]=v;
      offsets_rec(dim, d+1, cur, out);
    }
  }
}

// [[Rcpp::export]]
Rcpp::IntegerVector fh_meanshift_segmenter_cpp(
    NumericVector img_d,
    int nrow, int ncol, int nb,
    double fh_k = 0.5,
    int fh_min_size = 20,
    bool eight = true,
    int ms_dim = 0,                 // NEW: number of spectral dims to use (<=3 recommended)
    double ms_ranger = 0.2,
    double ms_hs = 10.0,
    int ms_max_iter = 10,
    double ms_eps = 1e-3,
    double mode_merge = 0.6,
    int final_min_size = 50,
    bool keep_invalid_na = true
) {
  if(nrow<=0 || ncol<=0 || nb<=0) stop("nrow/ncol/nb must be > 0");
  const int npix = nrow*ncol;
  if((int)img_d.size() != npix*nb) stop("img length must be nrow*ncol*nb");
  if(fh_k <= 0) stop("fh_k must be > 0");
  if(ms_ranger <= 0) stop("ms_ranger must be > 0");
  if(ms_hs <= 0) stop("ms_hs must be > 0");
  if(fh_min_size < 1) fh_min_size = 1;
  if(final_min_size < 1) final_min_size = 1;
  if(ms_max_iter < 1) ms_max_iter = 1;

  // Choose spectral feature dimension: keep dim<=5 => ms_dim + 2 <= 5 => ms_dim <= 3
  if(ms_dim <= 0) ms_dim = std::min(3, nb);
  ms_dim = std::max(1, std::min(ms_dim, std::min(3, nb)));

  const int dim = ms_dim + 2; // spectral + (row,col)
  // dim is now guaranteed in [3..5]

  // Convert to float, band-major
  std::vector<float> img((size_t)npix*nb);
  for(int i=0;i<npix*nb;++i) img[i] = (float)img_d[i];

  // Valid mask
  std::vector<unsigned char> valid(npix, 1);
  for(int p=0;p<npix;++p){
    for(int k=0;k<nb;++k){
      double v = img_d[p + npix*k];
      if(NumericVector::is_na(v) || !std::isfinite(v)){
        valid[p]=0; break;
      }
    }
  }

  auto idx = [ncol](int r,int c){ return r*ncol + c; };

  // -----------------------------
  // A) FH superpixels
  // -----------------------------
  std::vector<Edge> edges;
  edges.reserve((size_t)npix * 4);

  for(int r=0;r<nrow;++r){
    for(int c=0;c<ncol;++c){
      int a = idx(r,c);
      if(!valid[a]) continue;

      if(c+1<ncol){
        int b = idx(r,c+1);
        if(valid[b]) edges.push_back({a,b, pixDist(img.data(), npix, nb, a, b)});
      }
      if(r+1<nrow){
        int b = idx(r+1,c);
        if(valid[b]) edges.push_back({a,b, pixDist(img.data(), npix, nb, a, b)});
      }
      if(eight){
        if(r+1<nrow && c+1<ncol){
          int b = idx(r+1,c+1);
          if(valid[b]) edges.push_back({a,b, pixDist(img.data(), npix, nb, a, b)});
        }
        if(r+1<nrow && c-1>=0){
          int b = idx(r+1,c-1);
          if(valid[b]) edges.push_back({a,b, pixDist(img.data(), npix, nb, a, b)});
        }
      }
    }
  }

  std::sort(edges.begin(), edges.end(),
            [](const Edge& e1, const Edge& e2){ return e1.w < e2.w; });

  DSUFH dsu(npix);
  for(int p=0;p<npix;++p) if(!valid[p]) dsu.sz[p]=0;

  for(const auto& e: edges){
    int a = dsu.find(e.a);
    int b = dsu.find(e.b);
    if(a==b) continue;
    if(dsu.sz[a]==0 || dsu.sz[b]==0) continue;

    float ta = tau(dsu.internal[a], dsu.sz[a], (float)fh_k);
    float tb = tau(dsu.internal[b], dsu.sz[b], (float)fh_k);
    if(e.w <= std::min(ta,tb)){
      dsu.unite(a,b,e.w);
    }
  }

  if(fh_min_size > 1){
    for(const auto& e: edges){
      int a = dsu.find(e.a);
      int b = dsu.find(e.b);
      if(a==b) continue;
      if(dsu.sz[a]==0 || dsu.sz[b]==0) continue;
      if(dsu.sz[a] < fh_min_size || dsu.sz[b] < fh_min_size){
        dsu.unite(a,b,e.w);
      }
    }
  }

  // Remap FH components to 0..R-1
  std::unordered_map<int,int> comp_map;
  comp_map.reserve((size_t)npix/4);

  std::vector<int> fh_label(npix, -1);
  int R = 0;
  for(int p=0;p<npix;++p){
    if(!valid[p]) continue;
    int rp = dsu.find(p);
    auto it = comp_map.find(rp);
    if(it==comp_map.end()){
      comp_map.emplace(rp, R);
      fh_label[p] = R;
      R++;
    } else {
      fh_label[p] = it->second;
    }
  }

  if(R==0){
    IntegerVector out(npix, keep_invalid_na ? NA_INTEGER : 0);
    return out;
  }

  // -----------------------------
  // B) Region stats: mean spectra (first ms_dim bands) + centroid
  // -----------------------------
  std::vector<double> sumSpec((size_t)R * ms_dim, 0.0);
  std::vector<double> sumRow(R, 0.0), sumCol(R, 0.0);
  std::vector<int>    rSize(R, 0);

  for(int r=0;r<nrow;++r){
    for(int c=0;c<ncol;++c){
      int p = idx(r,c);
      if(!valid[p]) continue;
      int rid = fh_label[p];
      rSize[rid] += 1;
      sumRow[rid] += (double)r;
      sumCol[rid] += (double)c;
      for(int k=0;k<ms_dim;++k){
        sumSpec[(size_t)rid*ms_dim + k] += (double)img[p + npix*k];
      }
    }
  }

  // Build normalized feature z = [meanSpec/ms_ranger, row/ms_hs, col/ms_hs]
  std::vector<float> z((size_t)R * dim);

  const float invRng = (float)(1.0 / ms_ranger);
  const float invHs  = (float)(1.0 / ms_hs);

  for(int i=0;i<R;++i){
    const double inv = 1.0 / std::max(1, rSize[i]);
    for(int d=0; d<ms_dim; ++d){
      z[(size_t)i*dim + d] = (float)(sumSpec[(size_t)i*ms_dim + d] * inv) * invRng;
    }
    z[(size_t)i*dim + ms_dim + 0] = (float)(sumRow[i] * inv) * invHs;
    z[(size_t)i*dim + ms_dim + 1] = (float)(sumCol[i] * inv) * invHs;
  }

  // -----------------------------
  // C) MeanShift in region feature space (grid hash)
  // -----------------------------
  std::unordered_map<Key, std::vector<int>, KeyHash> grid;
  grid.reserve((size_t)R*2);

  for(int i=0;i<R;++i){
    Key k = makeKey(z, i, dim);
    grid[k].push_back(i);
  }

  std::vector< std::vector<int> > neigh;
  std::vector<int> cur(dim, 0);
  offsets_rec(dim, 0, cur, neigh);

  std::vector<float> y = z;
  std::vector<float> ynew((size_t)R * dim);

  const float eps2 = (float)(ms_eps * ms_eps);

  auto kgauss = [](float d2){ return std::exp(-0.5f * d2); };

  for(int i=0;i<R;++i){
    for(int it=0; it<ms_max_iter; ++it){
      Key qk; qk.dim=dim;
      qk.x0 = (int)std::floor(y[(size_t)i*dim + 0]);
      qk.x1 = (int)std::floor(y[(size_t)i*dim + 1]);
      qk.x2 = (dim>2)? (int)std::floor(y[(size_t)i*dim + 2]) : 0;
      qk.x3 = (dim>3)? (int)std::floor(y[(size_t)i*dim + 3]) : 0;
      qk.x4 = (dim>4)? (int)std::floor(y[(size_t)i*dim + 4]) : 0;

      float wsum_f = 0.0f;
      std::vector<float> acc(dim, 0.0f);

      for(const auto& off : neigh){
        Key kk = qk;
        kk.x0 += off[0]; kk.x1 += off[1];
        if(dim>2) kk.x2 += off[2];
        if(dim>3) kk.x3 += off[3];
        if(dim>4) kk.x4 += off[4];

        auto itg = grid.find(kk);
        if(itg == grid.end()) continue;

        for(int j : itg->second){
          float d2 = 0.0f;
          for(int d=0; d<dim; ++d){
            float diff = z[(size_t)j*dim + d] - y[(size_t)i*dim + d];
            d2 += diff*diff;
          }
          if(d2 > 9.0f) continue; // 3-sigma cut (in normalized feature units)
          float w = kgauss(d2) * (float)rSize[j]; // size-weighted
          wsum_f += w;
          for(int d=0; d<dim; ++d) acc[d] += w * z[(size_t)j*dim + d];
        }
      }

      if(wsum_f <= 0.0f) break;

      float shift2 = 0.0f;
      for(int d=0; d<dim; ++d){
        float m = acc[d] / wsum_f;
        float diff = m - y[(size_t)i*dim + d];
        shift2 += diff*diff;
        ynew[(size_t)i*dim + d] = m;
      }

      for(int d=0; d<dim; ++d) y[(size_t)i*dim + d] = ynew[(size_t)i*dim + d];
      if(shift2 <= eps2) break;
    }
  }

  // -----------------------------
  // D) Mode clustering (merge nearby converged points)
  // -----------------------------
  const float mode2 = (float)(mode_merge * mode_merge);

  std::vector<float> ym((size_t)R * dim);
  for(int i=0;i<R;++i){
    for(int d=0; d<dim; ++d) ym[(size_t)i*dim + d] = y[(size_t)i*dim + d] / (float)mode_merge;
  }

  std::unordered_map<Key, std::vector<int>, KeyHash> modeGrid;
  modeGrid.reserve((size_t)R*2);
  for(int i=0;i<R;++i){
    Key k = makeKey(ym, i, dim);
    modeGrid[k].push_back(i);
  }

  DSU modeDSU(R);
  for(int i=0;i<R;++i) modeDSU.sz[i] = rSize[i];

  std::vector< std::vector<int> > neigh2;
  std::vector<int> cur2(dim, 0);
  offsets_rec(dim, 0, cur2, neigh2);

  for(int i=0;i<R;++i){
    Key qk = makeKey(ym, i, dim);
    for(const auto& off : neigh2){
      Key kk = qk;
      kk.x0 += off[0]; kk.x1 += off[1];
      if(dim>2) kk.x2 += off[2];
      if(dim>3) kk.x3 += off[3];
      if(dim>4) kk.x4 += off[4];

      auto itg = modeGrid.find(kk);
      if(itg == modeGrid.end()) continue;

      for(int j : itg->second){
        if(i==j) continue;
        float d2=0.0f;
        for(int d=0; d<dim; ++d){
          float diff = y[(size_t)i*dim + d] - y[(size_t)j*dim + d];
          d2 += diff*diff;
        }
        if(d2 <= mode2) modeDSU.unite(i,j);
      }
    }
  }

  // -----------------------------
  // E) Region adjacency from FH labels (unique edges)
  // -----------------------------
  std::unordered_set<std::uint64_t> adjSet;
  adjSet.reserve((size_t)R*4);

  for(int r=0;r<nrow;++r){
    for(int c=0;c<ncol;++c){
      int p = idx(r,c);
      if(!valid[p]) continue;
      int a = fh_label[p];

      if(c+1<ncol){
        int q = idx(r,c+1);
        if(valid[q]){
          int b = fh_label[q];
          if(a!=b) adjSet.insert(packEdge(a,b));
        }
      }
      if(r+1<nrow){
        int q = idx(r+1,c);
        if(valid[q]){
          int b = fh_label[q];
          if(a!=b) adjSet.insert(packEdge(a,b));
        }
      }
      if(eight){
        if(r+1<nrow && c+1<ncol){
          int q = idx(r+1,c+1);
          if(valid[q]){
            int b = fh_label[q];
            if(a!=b) adjSet.insert(packEdge(a,b));
          }
        }
        if(r+1<nrow && c-1>=0){
          int q = idx(r+1,c-1);
          if(valid[q]){
            int b = fh_label[q];
            if(a!=b) adjSet.insert(packEdge(a,b));
          }
        }
      }
    }
  }

  std::vector<std::pair<int,int>> adj;
  adj.reserve(adjSet.size());
  for(auto key : adjSet){
    int a = (int)(key >> 32);
    int b = (int)(key & 0xFFFFFFFFu);
    adj.emplace_back(a,b);
  }

  // Map FH region -> mode root
  std::vector<int> reg2mode(R);
  for(int i=0;i<R;++i) reg2mode[i] = modeDSU.find(i);

  // Mode sizes (pixel counts)
  std::unordered_map<int,int> modeSize;
  modeSize.reserve((size_t)R);
  for(int i=0;i<R;++i) modeSize[reg2mode[i]] += rSize[i];

  // Compress mode roots to 0..nModes-1
  std::unordered_map<int,int> modeMap;
  modeMap.reserve(modeSize.size());
  int nModes = 0;
  for(const auto& kv : modeSize){
    modeMap[kv.first] = nModes++;
  }

  DSU finalDSU(nModes);
  for(const auto& kv : modeSize){
    int mid = modeMap[kv.first];
    finalDSU.sz[mid] = kv.second;
  }

  // Build edges between modes using region-feature distance proxy
  struct MEdge { int a,b; float d2; };
  std::vector<MEdge> medges;
  medges.reserve(adj.size());

  for(const auto& e : adj){
    int ra = e.first, rb = e.second;
    int ma = modeMap[ reg2mode[ra] ];
    int mb = modeMap[ reg2mode[rb] ];
    if(ma==mb) continue;

    float d2 = 0.0f;
    for(int d=0; d<dim; ++d){
      float diff = y[(size_t)ra*dim + d] - y[(size_t)rb*dim + d];
      d2 += diff*diff;
    }
    medges.push_back({ma,mb,d2});
  }

  std::sort(medges.begin(), medges.end(),
            [](const MEdge& x, const MEdge& y){ return x.d2 < y.d2; });

  if(final_min_size > 1){
    for(const auto& e : medges){
      int a = finalDSU.find(e.a);
      int b = finalDSU.find(e.b);
      if(a==b) continue;
      if(finalDSU.sz[a] < final_min_size || finalDSU.sz[b] < final_min_size){
        finalDSU.unite(a,b);
      }
    }
  }

  // Final relabel of mode DSU roots to 1..K
  std::unordered_map<int,int> finalRemap;
  finalRemap.reserve((size_t)nModes);

  std::vector<int> modeToFinal(nModes, 0);
  int nextLab = 1;
  for(int m=0; m<nModes; ++m){
    int rm = finalDSU.find(m);
    auto it = finalRemap.find(rm);
    if(it==finalRemap.end()){
      finalRemap[rm] = nextLab;
      modeToFinal[m] = nextLab++;
    } else {
      modeToFinal[m] = it->second;
    }
  }

  IntegerVector out(npix);
  for(int p=0;p<npix;++p){
    if(!valid[p]){
      out[p] = keep_invalid_na ? NA_INTEGER : 0;
      continue;
    }
    int rid = fh_label[p];
    int mid = modeMap[ reg2mode[rid] ];
    out[p] = modeToFinal[ finalDSU.find(mid) ];
  }

  return out;
}
