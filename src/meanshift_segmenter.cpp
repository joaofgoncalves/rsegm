#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>

using namespace Rcpp;

namespace{

  // ---------------------------
  // Union-Find (DSU)
  // ---------------------------
  struct DSU {
    std::vector<int> p, r, sz;
    DSU(int n): p(n), r(n,0), sz(n,1) { for(int i=0;i<n;++i) p[i]=i; }
    inline int find(int x){
      while(p[x]!=x){ p[x]=p[p[x]]; x=p[x]; }
      return x;
    }
    inline void unite(int a,int b){
      a=find(a); b=find(b);
      if(a==b) return;
      if(r[a]<r[b]) std::swap(a,b);
      p[b]=a; sz[a]+=sz[b];
      if(r[a]==r[b]) r[a]++;
    }
  };

  // squared distance between vectors of length nb stored in band-major layout
  inline float dist2_bandmajor(const float* img, int npix, int nb, int p, const float* y){
    float acc=0.0f;
    for(int k=0;k<nb;++k){
      float d = img[p + npix*k] - y[k];
      acc += d*d;
    }
    return acc;
  }

  inline float dist2_filtered(const float* filt, int npix, int nb, int p, int q){
    float acc=0.0f;
    for(int k=0;k<nb;++k){
      float d = filt[p + npix*k] - filt[q + npix*k];
      acc += d*d;
    }
    return acc;
  }


}


// [[Rcpp::export]]
Rcpp::List meanshift_segmenter_cpp(NumericVector img_d,
                           int nrow, int ncol, int nb,
                           int spatialr = 5,
                           double ranger = 1.0,
                           int max_iter = 10,
                           double eps = 1e-3,
                           double merge_thr = NA_REAL,
                           int min_size = 30,
                           bool eight = true) {
  const int npix = nrow*ncol;
  if((int)img_d.size() != npix*nb) stop("img length must be nrow*ncol*nb");
  if(spatialr < 1) spatialr = 1;
  if(ranger <= 0) stop("ranger must be > 0");
  if(max_iter < 1) max_iter = 1;
  if(min_size < 1) min_size = 1;

  // Convert to float for speed
  std::vector<float> img(npix*nb);
  for(int i=0;i<npix*nb;++i) img[i] = (float)img_d[i];

  // Valid mask: pixel valid if all bands finite & not NA
  std::vector<unsigned char> valid(npix, 1);
  for(int p=0;p<npix;++p){
    for(int k=0;k<nb;++k){
      double v = img_d[p + npix*k];
      if(NumericVector::is_na(v) || !std::isfinite(v)){
        valid[p]=0; break;
      }
    }
  }

  // Precompute spatial kernel weights for offsets in [-spatialr, spatialr]
  // Gaussian spatial kernel with sigma = spatialr/2 (good default)
  const float sig_s = std::max(1.0f, spatialr/2.0f);
  const float inv2sig_s2 = 1.0f / (2.0f*sig_s*sig_s);

  std::vector<int> dR;
  std::vector<int> dC;
  std::vector<float> wS;
  dR.reserve((2*spatialr+1)*(2*spatialr+1));
  dC.reserve((2*spatialr+1)*(2*spatialr+1));
  wS.reserve((2*spatialr+1)*(2*spatialr+1));

  for(int dr=-spatialr; dr<=spatialr; ++dr){
    for(int dc=-spatialr; dc<=spatialr; ++dc){
      if(dr==0 && dc==0) continue;
      if(!eight && std::abs(dr)+std::abs(dc)==2) continue;
      float ds2 = (float)(dr*dr + dc*dc);
      float ws = std::exp(-ds2 * inv2sig_s2);
      dR.push_back(dr); dC.push_back(dc); wS.push_back(ws);
    }
  }

  const float hr = (float)ranger;
  const float inv2hr2 = 1.0f / (2.0f*hr*hr);
  const float eps2 = (float)(eps*eps);

  // Output filtered image (band-major)
  std::vector<float> filt(npix*nb, 0.0f);

  // Mean-shift filtering per pixel (fixed samples: img)
  // For stability + speed: start at original value
  std::vector<float> y(nb), ynew(nb);

  auto idx = [ncol](int r,int c){ return r*ncol + c; };

  for(int r=0; r<nrow; ++r){
    for(int c=0; c<ncol; ++c){
      int p = idx(r,c);
      if(!valid[p]){
        // keep NA as 0 in filt; label will become 0 later
        continue;
      }

      // init y = x_p
      for(int k=0;k<nb;++k) y[k] = img[p + npix*k];

      for(int it=0; it<max_iter; ++it){
        // weighted mean
        float wsum = 1.0f;
        for(int k=0;k<nb;++k) ynew[k] = y[k]; // include self with weight 1

        // neighborhood accumulation
        for(size_t t=0; t<dR.size(); ++t){
          int rr = r + dR[t];
          int cc = c + dC[t];
          if(rr<0 || cc<0 || rr>=nrow || cc>=ncol) continue;
          int q = idx(rr,cc);
          if(!valid[q]) continue;

          // range weight
          float d2 = dist2_bandmajor(img.data(), npix, nb, q, y.data());
          float wr = std::exp(-d2 * inv2hr2);
          float w = wS[t] * wr;
          if(w < 1e-6f) continue;

          wsum += w;
          for(int k=0;k<nb;++k){
            ynew[k] += w * img[q + npix*k];
          }
        }

        float inv = 1.0f / wsum;
        float shift2 = 0.0f;
        for(int k=0;k<nb;++k){
          float m = ynew[k] * inv;
          float d = m - y[k];
          shift2 += d*d;
          y[k] = m;
        }
        if(shift2 <= eps2) break;
      }

      // store filtered
      for(int k=0;k<nb;++k) filt[p + npix*k] = y[k];
    }
  }

  // Default merge threshold: 0.5 * ranger (common practical choice)
  double mthr = Rcpp::NumericVector::is_na(merge_thr) ? (0.5 * ranger) : merge_thr;
  if(mthr <= 0) mthr = 0.5 * ranger;
  const float mthr2 = (float)(mthr*mthr);

  // Step B: connected-component clustering based on filtered distance
  DSU dsu(npix);
  for(int p=0;p<npix;++p) if(!valid[p]) dsu.sz[p]=0;

  for(int r=0; r<nrow; ++r){
    for(int c=0; c<ncol; ++c){
      int p = idx(r,c);
      if(!valid[p]) continue;

      // right, down (+ diagonals if eight)
      if(c+1<ncol){
        int q=idx(r,c+1);
        if(valid[q] && dist2_filtered(filt.data(), npix, nb, p, q) <= mthr2) dsu.unite(p,q);
      }
      if(r+1<nrow){
        int q=idx(r+1,c);
        if(valid[q] && dist2_filtered(filt.data(), npix, nb, p, q) <= mthr2) dsu.unite(p,q);
      }
      if(eight){
        if(r+1<nrow && c+1<ncol){
          int q=idx(r+1,c+1);
          if(valid[q] && dist2_filtered(filt.data(), npix, nb, p, q) <= mthr2) dsu.unite(p,q);
        }
        if(r+1<nrow && c-1>=0){
          int q=idx(r+1,c-1);
          if(valid[q] && dist2_filtered(filt.data(), npix, nb, p, q) <= mthr2) dsu.unite(p,q);
        }
      }
    }
  }

  // Step C: min_size enforcement by merging small components to most similar neighbor
  if(min_size > 1){
    // For each pixel, if its component is small, try to merge via best neighbor edge (local)
    for(int r=0; r<nrow; ++r){
      for(int c=0; c<ncol; ++c){
        int p = idx(r,c);
        if(!valid[p]) continue;
        int rp = dsu.find(p);
        if(dsu.sz[rp] >= min_size) continue;

        float best = 1e30f;
        int bestRoot = -1;

        // search local neighbors (same list as spatial offsets, plus self not needed)
        for(size_t t=0; t<dR.size(); ++t){
          int rr = r + dR[t];
          int cc = c + dC[t];
          if(rr<0 || cc<0 || rr>=nrow || cc>=ncol) continue;
          int q = idx(rr,cc);
          if(!valid[q]) continue;
          int rq = dsu.find(q);
          if(rq == rp || dsu.sz[rq]==0) continue;

          float d2 = dist2_filtered(filt.data(), npix, nb, p, q);
          if(d2 < best){
            best = d2;
            bestRoot = rq;
          }
        }

        if(bestRoot >= 0){
          dsu.unite(rp, bestRoot);
        }
      }
    }
  }

  // Relabel components 1..K; invalid -> 0
  IntegerVector labels(npix, 0);
  std::unordered_map<int,int> remap;
  remap.reserve((size_t)npix/4);
  int nextLab=1;

  for(int p=0;p<npix;++p){
    if(!valid[p]) { labels[p]=0; continue; }
    int rp = dsu.find(p);
    auto it = remap.find(rp);
    if(it==remap.end()){
      remap.emplace(rp, nextLab);
      labels[p]=nextLab++;
    } else {
      labels[p]=it->second;
    }
  }

  // Return both labels and filtered image (optional for debugging/inspection)
  NumericVector filt_out(npix*nb);
  for(int i=0;i<npix*nb;++i) filt_out[i] = (double)filt[i];

  return List::create(
    _["labels"] = labels,
    _["filtered"] = filt_out,
    _["nrow"] = nrow,
    _["ncol"] = ncol,
    _["nb"] = nb
  );
}
