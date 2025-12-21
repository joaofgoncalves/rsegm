// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <cmath>
#include <cstdint>
#include <algorithm>

using namespace Rcpp;

namespace{

static inline std::uint64_t pair_key(int a, int b) {
  if (a > b) std::swap(a,b);
  return ( (std::uint64_t)(std::uint32_t)a << 32 ) | (std::uint64_t)(std::uint32_t)b;
}
static inline double safe_sqrt(double x){ return std::sqrt(x > 0.0 ? x : 0.0); }

// DSU for pixel->final region mapping (cheap, necessary)
struct DSU {
  std::vector<int> p, r;
  DSU(int n): p(n), r(n,0){ for(int i=0;i<n;++i) p[i]=i; }
  inline int find(int x){
    while(p[x]!=x){ p[x]=p[p[x]]; x=p[x]; }
    return x;
  }
  inline int unite(int a,int b){
    a=find(a); b=find(b);
    if(a==b) return a;
    if(r[a]<r[b]) std::swap(a,b);
    p[b]=a;
    if(r[a]==r[b]) r[a]++;
    return a;
  }
};

struct Region {
  int id;
  bool alive;
  int n;
  int minr,maxr,minc,maxc;
  double perim;
  int version;
  std::vector<double> sum, sumsq;
  std::unordered_set<int> neigh;

  Region() = delete;
  Region(int id_, int nb, int r, int c):
    id(id_), alive(true), n(1),
    minr(r),maxr(r),minc(c),maxc(c),
    perim(4.0), version(0),
    sum(nb,0.0), sumsq(nb,0.0) {}
};

struct Cand {
  int a,b;
  double ff;
  int va,vb;
  bool operator>(Cand const& o) const { return ff > o.ff; }
};

static inline double spectral_increase(const Region& A, const Region& B,
                                       const std::vector<double>& bw){
  const int nb = (int)A.sum.size();
  const int n1=A.n, n2=B.n, n3=n1+n2;
  if(n3<=0) return 0.0;
  double h=0.0;
  for(int k=0;k<nb;++k){
    const double m1=A.sum[k]/n1, m2=B.sum[k]/n2, m3=(A.sum[k]+B.sum[k])/n3;
    const double v1=(A.sumsq[k]/n1)-m1*m1;
    const double v2=(B.sumsq[k]/n2)-m2*m2;
    const double v3=((A.sumsq[k]+B.sumsq[k])/n3)-m3*m3;
    const double sd1=safe_sqrt(v1), sd2=safe_sqrt(v2), sd3=safe_sqrt(v3);
    h += bw[k]*((double)n3*sd3 - (double)n1*sd1 - (double)n2*sd2);
  }
  return h;
}

static inline double shape_increase(const Region& A, const Region& B,
                                    double perim3, double bboxp3,
                                    double w_compact){
  const int n1=A.n, n2=B.n, n3=n1+n2;
  if(n3<=0) return 0.0;
  const double l1=A.perim, l2=B.perim, l3=perim3;
  const double b1=2.0*((A.maxc-A.minc+1)+(A.maxr-A.minr+1));
  const double b2=2.0*((B.maxc-B.minc+1)+(B.maxr-B.minr+1));
  const double b3=bboxp3;

  const double cmp1=(n1>0)? l1/std::sqrt((double)n1):0.0;
  const double cmp2=(n2>0)? l2/std::sqrt((double)n2):0.0;
  const double cmp3=(n3>0)? l3/std::sqrt((double)n3):0.0;
  const double h_cmp=(double)n3*cmp3 - (double)n1*cmp1 - (double)n2*cmp2;

  const double sm1=(n1>0 && b1>0)? l1/b1:0.0;
  const double sm2=(n2>0 && b2>0)? l2/b2:0.0;
  const double sm3=(n3>0 && b3>0)? l3/b3:0.0;
  const double h_sm=(double)n3*sm3 - (double)n1*sm1 - (double)n2*sm2;

  return w_compact*h_cmp + (1.0-w_compact)*h_sm;
}

static inline double fusion_factor(const Region& A, const Region& B,
                                   const std::vector<double>& bw,
                                   double w_color, double w_compact,
                                   int sharedAB){
  const double perim3 = A.perim + B.perim - 2.0*(double)sharedAB;
  const int minr=std::min(A.minr,B.minr), maxr=std::max(A.maxr,B.maxr);
  const int minc=std::min(A.minc,B.minc), maxc=std::max(A.maxc,B.maxc);
  const double bboxp3 = 2.0*((maxc-minc+1)+(maxr-minr+1));
  const double hc = spectral_increase(A,B,bw);
  const double hs = shape_increase(A,B,perim3,bboxp3,w_compact);
  return w_color*hc + (1.0-w_color)*hs;
}

}

// [[Rcpp::export]]
Rcpp::IntegerVector baatz_fast_segmenter_cpp(const arma::mat& image,
                                              int nrows,
                                              int ncols,
                                              double scale_param,
                                              double color_weight = 0.9,
                                              double compactness_weight = 0.5,
                                              Nullable<NumericVector> band_weights = R_NilValue,
                                              bool verbose = false) {
  const int np = nrows*ncols;
  const int nb = image.n_cols;
  if((int)image.n_rows != np) stop("image must have nrows*ncols rows");

  std::vector<double> bw(nb,1.0);
  if(band_weights.isNotNull()){
    NumericVector bwr(band_weights);
    if((int)bwr.size()!=nb) stop("band_weights length must equal nbands");
    for(int k=0;k<nb;++k) bw[k]=(double)bwr[k];
  }

  auto idx = [ncols](int r,int c){ return r*ncols + c; };

  std::vector<Region> R; R.reserve(np);
  for(int r=0;r<nrows;++r){
    for(int c=0;c<ncols;++c){
      int id=idx(r,c);
      R.emplace_back(id, nb, r, c);
      Region& reg = R.back();
      for(int k=0;k<nb;++k){
        double v = image(id,k);
        reg.sum[k]=v;
        reg.sumsq[k]=v*v;
      }
    }
  }

  DSU dsu(np);

  std::unordered_map<std::uint64_t,int> shared;
  shared.reserve((size_t)np*2);

  for(int r=0;r<nrows;++r){
    for(int c=0;c<ncols;++c){
      int a=idx(r,c);
      if(c+1<ncols){
        int b=idx(r,c+1);
        R[a].neigh.insert(b);
        R[b].neigh.insert(a);
        shared.emplace(pair_key(a,b), 1);
      }
      if(r+1<nrows){
        int b=idx(r+1,c);
        R[a].neigh.insert(b);
        R[b].neigh.insert(a);
        shared.emplace(pair_key(a,b), 1);
      }
    }
  }

  auto get_shared = [&](int a,int b)->int{
    auto it=shared.find(pair_key(a,b));
    return (it==shared.end())?0:it->second;
  };
  auto erase_shared = [&](int a,int b){ shared.erase(pair_key(a,b)); };
  auto add_shared = [&](int a,int b,int addlen){
    if(a==b) return;
    std::uint64_t k=pair_key(a,b);
    auto it=shared.find(k);
    if(it==shared.end()) shared.emplace(k, addlen);
    else it->second += addlen;
  };

  const double thr = scale_param*scale_param;

  std::priority_queue<Cand, std::vector<Cand>, std::greater<Cand>> pq;

  for(const auto& kv: shared){
    std::uint64_t key=kv.first;
    int a=(int)(key>>32);
    int b=(int)(key & 0xFFFFFFFFu);
    int sh=kv.second;
    double ff=fusion_factor(R[a],R[b],bw,color_weight,compactness_weight,sh);
    pq.push({a,b,ff,R[a].version,R[b].version});
  }

  int alive=np;
  int merges=0;

  auto merge_into = [&](int A,int B){
    // DSU: attach B under A (note: region roots are pixel ids; OK for tile sizes)
    dsu.unite(A,B);

    Region& ra=R[A];
    Region& rb=R[B];

    const int shAB=get_shared(A,B);

    ra.perim = ra.perim + rb.perim - 2.0*(double)shAB;

    ra.minr=std::min(ra.minr,rb.minr);
    ra.maxr=std::max(ra.maxr,rb.maxr);
    ra.minc=std::min(ra.minc,rb.minc);
    ra.maxc=std::max(ra.maxc,rb.maxc);

    ra.n += rb.n;
    for(int k=0;k<nb;++k){
      ra.sum[k]+=rb.sum[k];
      ra.sumsq[k]+=rb.sumsq[k];
    }

    for(int nbh: rb.neigh){
      if(!R[nbh].alive) continue;
      if(nbh==A) continue;

      const int shBN=get_shared(B,nbh);
      if(shBN>0){
        add_shared(A,nbh,shBN);
        erase_shared(B,nbh);
      }
      R[nbh].neigh.erase(B);
      R[nbh].neigh.insert(A);
      ra.neigh.insert(nbh);
    }

    erase_shared(A,B);
    ra.neigh.erase(B);

    rb.alive=false;
    rb.neigh.clear();
    rb.n=0;

    ra.version++;
  };

  while(!pq.empty()){
    Cand c=pq.top(); pq.pop();

    int a=c.a, b=c.b;
    if(a==b) continue;
    if(!R[a].alive || !R[b].alive) continue;
    if(R[a].version!=c.va || R[b].version!=c.vb) continue;
    if(R[a].neigh.find(b)==R[a].neigh.end()) continue;

    int sh=get_shared(a,b);
    if(sh<=0) continue;

    double ff=fusion_factor(R[a],R[b],bw,color_weight,compactness_weight,sh);
    if(ff!=c.ff){
      pq.push({a,b,ff,R[a].version,R[b].version});
      continue;
    }

    if(ff>=thr) break;

    // merge smaller into larger for speed
    if(R[a].n < R[b].n) std::swap(a,b);

    merge_into(a,b);
    alive--; merges++;

    Region& ra=R[a];
    for(int nbh: ra.neigh){
      if(!R[nbh].alive || nbh==a) continue;
      int sh2=get_shared(a,nbh);
      if(sh2<=0) continue;
      double ff2=fusion_factor(ra,R[nbh],bw,color_weight,compactness_weight,sh2);
      pq.push({a,nbh,ff2,ra.version,R[nbh].version});
    }

    if(verbose && (merges % 5000 == 0)){
      Rcpp::Rcout << "Merges: " << merges << " alive: " << alive << "\n";
    }
  }

  // Build final labels: pixel p belongs to DSU root find(p)
  // But DSU roots may reference dead nodes; we remap all unique roots to 1..K.
  IntegerVector lab(np);
  std::unordered_map<int,int> remap;
  remap.reserve((size_t)np/4);

  int next=1;
  for(int p=0;p<np;++p){
    int r = dsu.find(p);
    auto it=remap.find(r);
    if(it==remap.end()){
      remap.emplace(r,next);
      lab[p]=next;
      next++;
    } else {
      lab[p]=it->second;
    }
  }
  return lab;
}
