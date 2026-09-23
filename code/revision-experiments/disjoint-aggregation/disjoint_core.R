# Fixed disjoint blocks: no random or outcome-dependent subset selection.
fixed_blocks <- function(n, m, K=10L) {
  stopifnot(m>=1L, K>=1L, m*K<=n)
  lapply(seq_len(K), function(k) seq.int((k-1L)*m+1L,k*m))
}

disjoint_aggregate <- function(D,y,subsets,K_grid=c(1L,5L,10L)) {
  stopifnot(!anyDuplicated(unlist(subsets)),max(K_grid)<=length(subsets))
  fits <- lapply(subsets,function(ids) ols_influence(D,y,ids))
  out <- lapply(K_grid,function(K) {
    f <- fits[seq_len(K)]
    beta <- Reduce(`+`,lapply(f,`[[`,"beta"))/K
    # This is exactly the conditional covariance of the specified linearized
    # independent Gaussian multiplier, not the exact finite-sample sampling variance.
    V <- Reduce(`+`,lapply(f,function(a) crossprod(a$contribution)))/K^2
    se <- sqrt(diag(V))
    list(K=K,beta=beta,covariance=V,se=se,
      lower=beta-qnorm(.975)*se,upper=beta+qnorm(.975)*se,
      status=if(all(vapply(f,function(a)a$status=="ok",logical(1))))"ok"else"subset_fit_failure",
      min_rcond=min(vapply(f,`[[`,numeric(1),"rcond")))
  })
  list(aggregates=out,blocks=fits)
}
