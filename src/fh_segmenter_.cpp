#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_map>

using namespace Rcpp;

struct Edge {
  int a;
  int b;
  float w;
};

struct DSU {
  std::vector<int> parent;
  std::vector<int> rnk;
  std::vector<int> sz;
  std::vector<float> internal;

  DSU(int n): parent(n), rnk(n,0), sz(n,1), internal(n,0.0f) {
    for(int i=0;i<n;++i) parent[i]=i;
  }

  inline int find(int x){
    while(parent[x]!=x){
      parent[x]=parent[parent[x]];
      x=parent[x];
    }
    return x;
  }

  inline void unite(int a,int b,float w){
    a=find(a); b=find(b);
    if(a==b) return;
    if(rnk[a]<rnk[b]) std::swap(a,b);
    parent[b]=a;
    sz[a]+=sz[b];
    internal[a]=std::max(std::max(internal[a], internal[b]), w);
    if(rnk[a]==rnk[b]) rnk[a]++;
  }
};

inline float tau(float internal, int size, float k){
  return internal + k / (float)size;
}

inline float pixDist_noNA(const double* img, int p, int q, int nb, int npix){
  // assumes both pixels are valid (no NA in any band)
  double acc = 0.0;
  for(int k=0;k<nb;++k){
    double d = img[p + npix*k] - img[q + npix*k];
    acc += d*d;
  }
  return (float)std::sqrt(acc);
}

// [[Rcpp::export]]
IntegerVector segment_fh_cpp(NumericVector img,
                                int nrow, int ncol, int nb,
                                double k = 0.8,
                                int min_size = 30,
                                bool eight = true) {
  const int npix = nrow * ncol;
  if((int)img.size() != npix * nb) stop("img length must be nrow*ncol*nb");
  if(k <= 0) stop("k must be > 0");
  if(min_size < 1) min_size = 1;

  const float kf = (float)k;
  const double* pimg = REAL(img);

  // Build valid mask: pixel is valid iff all bands are finite and not NA
  std::vector<unsigned char> valid(npix, 1);
  for(int p=0;p<npix;++p){
    for(int b=0;b<nb;++b){
      double v = pimg[p + npix*b];
      if(Rcpp::NumericVector::is_na(v) || !std::isfinite(v)){
        valid[p] = 0;
        break;
      }
    }
  }

  // Build edges only between valid pixels
  std::vector<Edge> edges;
  edges.reserve((size_t)npix * 4);

  auto id = [ncol](int r,int c){ return r*ncol + c; };

  for(int r=0;r<nrow;++r){
    for(int c=0;c<ncol;++c){
      int a = id(r,c);
      if(!valid[a]) continue;

      if(c+1 < ncol){
        int b = id(r,c+1);
        if(valid[b]) edges.push_back({a,b, pixDist_noNA(pimg,a,b,nb,npix)});
      }
      if(r+1 < nrow){
        int b = id(r+1,c);
        if(valid[b]) edges.push_back({a,b, pixDist_noNA(pimg,a,b,nb,npix)});
      }
      if(eight){
        if(r+1<nrow && c+1<ncol){
          int b = id(r+1,c+1);
          if(valid[b]) edges.push_back({a,b, pixDist_noNA(pimg,a,b,nb,npix)});
        }
        if(r+1<nrow && c-1>=0){
          int b = id(r+1,c-1);
          if(valid[b]) edges.push_back({a,b, pixDist_noNA(pimg,a,b,nb,npix)});
        }
      }
    }
  }

  std::sort(edges.begin(), edges.end(),
            [](const Edge& e1, const Edge& e2){ return e1.w < e2.w; });

  DSU dsu(npix);

  // Set invalid pixels to size 0 so they never influence min_size logic
  for(int p=0;p<npix;++p){
    if(!valid[p]) dsu.sz[p] = 0;
  }

  // Pass 1: FH adaptive merging
  for(const auto& e : edges){
    int a = dsu.find(e.a);
    int b = dsu.find(e.b);
    if(a==b) continue;

    // both components are valid if their size>0
    if(dsu.sz[a]==0 || dsu.sz[b]==0) continue;

    float ta = tau(dsu.internal[a], dsu.sz[a], kf);
    float tb = tau(dsu.internal[b], dsu.sz[b], kf);

    if(e.w <= std::min(ta, tb)){
      dsu.unite(a, b, e.w);
    }
  }

  // Pass 2: enforce minimum size on valid components
  if(min_size > 1){
    for(const auto& e : edges){
      int a = dsu.find(e.a);
      int b = dsu.find(e.b);
      if(a==b) continue;

      if(dsu.sz[a]==0 || dsu.sz[b]==0) continue;

      if(dsu.sz[a] < min_size || dsu.sz[b] < min_size){
        dsu.unite(a, b, e.w);
      }
    }
  }

  // Relabel valid pixels 1..K; invalid pixels keep label 0
  IntegerVector out(npix, 0);
  std::unordered_map<int,int> remap;
  remap.reserve((size_t)npix / 4);

  int nextLab=1;
  for(int p=0;p<npix;++p){
    if(!valid[p]) { out[p] = 0; continue; }
    int rp = dsu.find(p);
    auto it = remap.find(rp);
    if(it==remap.end()){
      remap.emplace(rp,nextLab);
      out[p]=nextLab++;
    } else {
      out[p]=it->second;
    }
  }

  return out;
}
