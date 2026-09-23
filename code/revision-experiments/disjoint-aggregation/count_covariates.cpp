#include <Rcpp.h>
#include <vector>
#include <cstdint>
using namespace Rcpp;

// The graph-generation and edgewise common-neighbor loops are copied from the
// anchor count_subgraphs.cpp. Weighted order-3/order-5 bootstrap moments are
// unnecessary for aggregation and are omitted. RNG order is unchanged.
// [[Rcpp::export]]
List count_covariates_cpp(NumericVector xi, double rho) {
  RNGScope scope;
  const int n=xi.size(), nw=(n+63)/64;
  std::vector<uint64_t> bits((size_t)n*nw,0);
  std::vector<int> eu,ev;
  eu.reserve((size_t)(rho*n*(n-1.0)/2*1.15+1000)); ev.reserve(eu.capacity());
  NumericVector degree(n),tri_twice(n);
  for(int i=0;i<n;++i)for(int j=i+1;j<n;++j) {
    if(unif_rand()<rho*(xi[i]+xi[j])) {
      eu.push_back(i);ev.push_back(j); degree[i]+=1;degree[j]+=1;
      bits[(size_t)i*nw+j/64]|=(uint64_t(1)<<(j%64));
      bits[(size_t)j*nw+i/64]|=(uint64_t(1)<<(i%64));
    }
  }
  for(size_t e=0;e<eu.size();++e) {
    const int i=eu[e],j=ev[e]; int common=0;
    for(int w=0;w<nw;++w)
      common+=__builtin_popcountll(bits[(size_t)i*nw+w]&bits[(size_t)j*nw+w]);
    tri_twice[i]+=common;tri_twice[j]+=common;
  }
  return List::create(_["degree"]=degree,_["triangles"]=tri_twice/2.0,
                      _["edges"]=(double)eu.size());
}
