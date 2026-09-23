#include <Rcpp.h>
#include <vector>
#include <cstdint>
#include <cmath>
using namespace Rcpp;

// Generate one graph and return everything needed for rooted two-star and
// rooted-triangle moments. Bitsets make triangle enumeration practical at
// n=6000 without constructing an n by n R matrix.
// [[Rcpp::export]]
List count_subgraphs_cpp(NumericVector xi, double rho, NumericMatrix weights) {
  RNGScope scope;
  const int n = xi.size();
  const int q = weights.ncol();
  const int nw = (n + 63) / 64;
  std::vector<uint64_t> bits((size_t)n * nw, 0);
  std::vector<int> eu, ev;
  const double expected_edges = rho * n * (n - 1.0) / 2.0;
  eu.reserve((size_t)(expected_edges * 1.15 + 1000));
  ev.reserve(eu.capacity());
  NumericVector degree(n);

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double p = rho * (xi[i] + xi[j]);
      if (unif_rand() < p) {
        eu.push_back(i); ev.push_back(j);
        degree[i] += 1.0; degree[j] += 1.0;
        bits[(size_t)i * nw + j / 64] |= (uint64_t(1) << (j % 64));
        bits[(size_t)j * nw + i / 64] |= (uint64_t(1) << (i % 64));
      }
    }
  }

  const size_t E = eu.size();
  std::vector<unsigned short> support((size_t)n * n, 0);
  NumericVector tri_twice(n), shared_edge_pairs(n);
  for (size_t e = 0; e < E; ++e) {
    const int i = eu[e], j = ev[e];
    int common = 0;
    for (int w = 0; w < nw; ++w)
      common += __builtin_popcountll(bits[(size_t)i*nw+w] & bits[(size_t)j*nw+w]);
    support[(size_t)i*n+j] = support[(size_t)j*n+i] = (unsigned short)common;
    tri_twice[i] += common; tri_twice[j] += common;
    const double pairs = common * (common - 1.0) / 2.0;
    shared_edge_pairs[i] += pairs; shared_edge_pairs[j] += pairs;
  }
  NumericVector triangles(n), bowtie_center(n), triangle_leaf(n);
  NumericMatrix ts_inclusion(n, q), tri_inclusion(n, q);
  NumericVector ts2_inclusion(n);
  for (int i = 0; i < n; ++i) {
    triangles[i] = tri_twice[i] / 2.0;
    bowtie_center[i] = triangles[i] * (triangles[i] - 1.0) / 2.0 - shared_edge_pairs[i];
    const double rooted4 = degree[i] >= 4 ?
      degree[i]*(degree[i]-1)*(degree[i]-2)*(degree[i]-3)/24.0 : 0.0;
    ts2_inclusion[i] = rooted4;
    const double rooted2 = degree[i] * (degree[i]-1.0) / 2.0;
    for (int c = 0; c < q; ++c) ts_inclusion(i,c) = weights(i,c) * rooted2;
  }

  // Peripheral contributions for the order-3 two-star moments and the
  // order-5 squared two-star leading term.
  for (size_t e = 0; e < E; ++e) {
    const int i = eu[e], j = ev[e];
    for (int c = 0; c < q; ++c) {
      ts_inclusion(j,c) += weights(i,c) * (degree[i] - 1.0);
      ts_inclusion(i,c) += weights(j,c) * (degree[j] - 1.0);
    }
    if (degree[i] >= 4)
      ts2_inclusion[j] += (degree[i]-1)*(degree[i]-2)*(degree[i]-3)/6.0;
    if (degree[j] >= 4)
      ts2_inclusion[i] += (degree[j]-1)*(degree[j]-2)*(degree[j]-3)/6.0;
  }

  // Enumerate each triangle once. Besides weighted order-3 inclusion sums,
  // accumulate the number of bowties in which each vertex is a leaf.
  for (size_t e = 0; e < E; ++e) {
    const int i = eu[e], j = ev[e]; // i < j
    for (int word = j / 64; word < nw; ++word) {
      uint64_t both = bits[(size_t)i*nw+word] & bits[(size_t)j*nw+word];
      if (word == j / 64) {
        const int shift = (j % 64) + 1;
        both &= shift == 64 ? uint64_t(0) : (~uint64_t(0) << shift);
      }
      while (both) {
        const int bit = __builtin_ctzll(both);
        const int k = word * 64 + bit;
        both &= both - 1;
        if (k >= n) continue;
        for (int c = 0; c < q; ++c) {
          const double h = (weights(i,c) + weights(j,c) + weights(k,c)) / 3.0;
          tri_inclusion(i,c) += h;
          tri_inclusion(j,c) += h;
          tri_inclusion(k,c) += h;
        }
        const double cij = support[(size_t)i*n+j];
        const double cik = support[(size_t)i*n+k];
        const double cjk = support[(size_t)j*n+k];
        const double vi = triangles[i] - cij - cik + 1.0;
        const double vj = triangles[j] - cij - cjk + 1.0;
        const double vk = triangles[k] - cik - cjk + 1.0;
        triangle_leaf[j] += vi; triangle_leaf[k] += vi;
        triangle_leaf[i] += vj; triangle_leaf[k] += vj;
        triangle_leaf[i] += vk; triangle_leaf[j] += vk;
      }
    }
  }
  NumericVector tri2_inclusion(n);
  for (int i = 0; i < n; ++i) tri2_inclusion[i] = bowtie_center[i] + triangle_leaf[i];

  return List::create(
    _["degree"] = degree,
    _["triangles"] = triangles,
    _["bowtie_center"] = bowtie_center,
    _["two_star_inclusion"] = ts_inclusion,
    _["triangle_inclusion"] = tri_inclusion,
    _["two_star_squared_inclusion"] = ts2_inclusion,
    _["triangle_squared_inclusion"] = tri2_inclusion,
    _["edges"] = (double)E
  );
}
