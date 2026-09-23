# Equal-weight aggregation of independently drawn, fixed-size regression subsets.
# Graph covariates are constructed once from the full graph. No graph resampling.

load_anchor <- function(anchor_dir, compile=TRUE) {
  env <- new.env(parent=globalenv())
  wanted <- c("truth", "clean_covariates", "safe_solve", "moment_fit", "validate_counts")
  for (expr in parse(file.path(anchor_dir, "coverage_experiment.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted) eval(expr, env)
  }
  stopifnot(all(vapply(wanted, exists, logical(1), envir=env, inherits=FALSE)))
  if (compile) Rcpp::sourceCpp(file.path(anchor_dir,"count_subgraphs.cpp"),env=env)
  env
}

ols_influence <- function(D, y, ids, tol=1e-12) {
  d <- D[ids,,drop=FALSE]; yy <- y[ids]; m <- length(ids); p <- ncol(D)
  M <- crossprod(d)/m; rc <- rcond(M)
  if (!is.finite(rc) || rc<tol) return(list(status="singular",rcond=rc,
    beta=rep(NA_real_,p),contribution=matrix(NA_real_,m,p)))
  beta <- drop(solve(M,crossprod(d,yy)/m))
  residual <- drop(yy-d%*%beta)
  # contribution_i = psi_i/m; HC0, matching the anchor's Gaussian multiplier CI.
  contribution <- t(solve(M,t(d*residual)))/m
  list(status="ok",rcond=rc,beta=beta,contribution=contribution)
}

aggregate_subsets <- function(D,y,subsets,K_grid=c(1L,5L,10L,20L)) {
  n <- nrow(D); p <- ncol(D); m <- length(subsets[[1]])
  stopifnot(length(y)==n, max(K_grid)<=length(subsets), all(K_grid>=1),
    all(vapply(subsets,length,integer(1))==m),
    all(vapply(subsets,function(x) !anyDuplicated(x) && all(x>=1 & x<=n),logical(1))))
  node_sum <- matrix(0,n,p); independent_var <- numeric(p); beta_sum <- numeric(p)
  multiplicity <- integer(n); min_rcond <- Inf; failures <- 0L; out <- list()
  for (k in seq_len(max(K_grid))) {
    ids <- subsets[[k]]; fit <- ols_influence(D,y,ids)
    min_rcond <- min(min_rcond,fit$rcond); failures <- failures+(fit$status!="ok")
    beta_sum <- beta_sum+fit$beta
    node_sum[ids,] <- node_sum[ids,,drop=FALSE]+fit$contribution
    independent_var <- independent_var+colSums(fit$contribution^2)
    multiplicity[ids] <- multiplicity[ids]+1L
    if (k %in% K_grid) {
      beta <- beta_sum/k
      # Identical node multipliers across subsets: cross terms retained exactly.
      shared_var <- colSums(node_sum^2)/k^2
      se <- sqrt(shared_var); naive_se <- sqrt(independent_var)/k
      out[[as.character(k)]] <- list(K=k,beta=beta,se=se,naive_se=naive_se,
        lower=beta-qnorm(.975)*se,upper=beta+qnorm(.975)*se,
        naive_lower=beta-qnorm(.975)*naive_se,naive_upper=beta+qnorm(.975)*naive_se,
        covariance=crossprod(node_sum)/k^2,
        overlap_variance=shared_var-independent_var/k^2,
        unique_nodes=sum(multiplicity>0),max_multiplicity=max(multiplicity),
        mean_pair_overlap=if(k>1)sum(choose(multiplicity,2))/choose(k,2) else NA_real_,
        status=if(failures)"subset_fit_failure" else "ok", failed_subsets=failures,
        min_rcond=min_rcond)
    }
  }
  out
}

summarize_aggregation <- function(x) {
  keys <- c("motif","delta","n","m","K","coefficient","method")
  groups <- split(x,interaction(x[keys],drop=TRUE,lex.order=TRUE))
  out <- lapply(groups,function(a) {
    ok <- a$status=="ok" & is.finite(a$estimate) & is.finite(a$lower) & is.finite(a$upper)
    b <- a[ok,,drop=FALSE]; N <- nrow(b); err <- b$estimate-b$target
    coverage <- if(N)mean(b$covered)else NA_real_
    data.frame(a[1,keys,drop=FALSE],attempted=nrow(a),valid=N,failed=nrow(a)-N,
      mean_estimate=if(N)mean(b$estimate)else NA_real_,
      bias=if(N)mean(err)else NA_real_,bias_mcse=if(N>1)sd(err)/sqrt(N)else NA_real_,
      empirical_sd=if(N>1)sd(b$estimate)else NA_real_,
      rmse=if(N)sqrt(mean(err^2))else NA_real_,
      mean_se=if(N)mean(b$se)else NA_real_,rms_se=if(N)sqrt(mean(b$se^2))else NA_real_,
      coverage=coverage,coverage_mcse=if(N)sqrt(coverage*(1-coverage)/N)else NA_real_,
      mean_width=if(N)mean(b$width)else NA_real_,
      width_mcse=if(N>1)sd(b$width)/sqrt(N)else NA_real_,
      mean_unique_nodes=mean(a$unique_nodes),
      mean_pair_overlap=if(all(is.na(a$mean_pair_overlap)))NA_real_ else mean(a$mean_pair_overlap),
      expected_pair_overlap=if(a$K[1]>1)a$m[1]^2/a$n[1]else NA_real_,
      mean_overlap_variance=mean(a$overlap_variance,na.rm=TRUE))
  })
  result <- do.call(rbind,out); rownames(result)<-NULL
  result[do.call(order,unname(as.list(result[keys]))),]
}

paired_changes <- function(x) {
  x <- x[x$method=="shared_node_multiplier",]
  join <- c("motif","delta","n","m","replicate","coefficient")
  base <- x[x$K==1,c(join,"estimate","covered","width")]
  names(base)[(length(join)+1):ncol(base)] <- c("base_estimate","base_covered","base_width")
  a <- merge(x[x$K>1,],base,by=join,all.x=TRUE)
  keys <- c("motif","delta","n","m","K","coefficient")
  groups <- split(a,interaction(a[keys],drop=TRUE,lex.order=TRUE))
  do.call(rbind,lapply(groups,function(b) {
    ok <- is.finite(b$estimate) & is.finite(b$base_estimate) & is.finite(b$width)
    b <- b[ok,]; d <- b$estimate-b$base_estimate
    dc <- as.numeric(b$covered)-as.numeric(b$base_covered); dw <- b$width-b$base_width
    data.frame(b[1,keys,drop=FALSE],paired_valid=nrow(b),
      estimate_change=mean(d),estimate_change_mcse=sd(d)/sqrt(length(d)),
      coverage_change=mean(dc),coverage_change_mcse=sd(dc)/sqrt(length(dc)),
      width_change=mean(dw),width_change_mcse=sd(dw)/sqrt(length(dw)))
  }))
}
