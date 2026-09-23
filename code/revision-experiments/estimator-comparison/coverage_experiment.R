suppressPackageStartupMessages(library(Rcpp))

script_arg <- grep("^--file=", commandArgs(), value=TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else "."
sourceCpp(file.path(script_dir, "count_subgraphs.cpp"))

truth <- c(intercept=1, network=20, x1=3, x2=2)
main_m_grid <- c(25L,40L,65L,100L,160L,250L,400L,650L,1000L)

clean_covariates <- function(x) {
  list(two_star=(x+.5)^2,
       triangle=x^2+7*x/6+1/3)
}

safe_solve <- function(M,b,tol=1e-12) {
  rc <- tryCatch(rcond(M),error=function(e) NA_real_)
  if(!is.finite(rc)||rc<tol) return(list(beta=rep(NA_real_,4),status="singular",rcond=rc))
  beta <- tryCatch(drop(solve(M,b)),error=function(e) rep(NA_real_,4))
  list(beta=beta,status=if(all(is.finite(beta))) "ok" else "solve_failure",rcond=rc)
}

moment_fit <- function(z,y,X,q=z^2,ids=seq_along(z)) {
  z <- z[ids]; y <- y[ids]; X <- X[ids,,drop=FALSE]; q <- q[ids]
  D <- cbind(1,z,X)
  Ai <- array(0,c(4,4,length(ids)))
  for(i in seq_along(ids)) {
    Ai[,,i] <- tcrossprod(D[i,])
    Ai[2,2,i] <- q[i]
  }
  bi <- D*y
  M <- apply(Ai,c(1,2),mean); b <- colMeans(bi)
  fit <- safe_solve(M,b)
  if(fit$status!="ok") return(c(fit,list(se=rep(NA_real_,4),lower=rep(NA_real_,4),upper=rep(NA_real_,4))))
  score <- bi
  for(i in seq_along(ids)) score[i,] <- bi[i,]-Ai[,,i]%*%fit$beta
  influence <- t(solve(M,t(score)))
  # Conditional Gaussian multiplier distribution is exactly normal.
  se <- sqrt(colSums(influence^2))/length(ids)
  crit <- qnorm(.975)
  c(fit,list(se=se,lower=fit$beta-crit*se,upper=fit$beta+crit*se))
}

full_structural_bootstrap <- function(motif,z,y,X,qcorr,inc,inc2,degree,rho_hat,
                                      B=999L,seed=1L) {
  n <- length(z); C2 <- choose(n-1,2); C4 <- choose(n-1,4)
  edges_per_copy <- if(motif=="two_star") 2 else 3
  divisor <- if(motif=="two_star") 3*C2*rho_hat^2 else C2*rho_hat^3
  zmean <- mean(z); yz <- mean(y*z); x1z <- mean(X[,1]*z); x2z <- mean(X[,2]*z)
  z2corr <- mean(qcorr); z2raw <- mean(z^2)
  gnet <- sweep(inc/divisor,2,c(zmean,yz,x1z,x2z),"-")
  divisor2 <- if(motif=="two_star") 5*C4*rho_hat^4 else 15*C4*rho_hat^6
  gz2 <- inc2/divisor2-z2corr
  rho_g <- degree/(n-1)-rho_hat

  ordinary <- cbind(y,X[,1],X[,2],X[,1]^2,X[,2]^2,X[,1]*X[,2],y*X[,1],y*X[,2])
  ordinary_mean <- colMeans(ordinary)
  ordinary_g <- sweep(ordinary,2,ordinary_mean,"-")
  set.seed(seed)
  W <- matrix(rnorm(n*B),n,B)
  net_delta <- crossprod(gnet,W)*(3/n)
  z2_delta <- drop(crossprod(gz2,W))*(5/n)
  ordinary_delta <- crossprod(ordinary_g,W)/n
  rho_b <- rho_hat+drop(crossprod(rho_g,W))*(2/n)
  ratio1 <- (rho_hat/rho_b)^edges_per_copy
  ratio2 <- (rho_hat/rho_b)^(2*edges_per_copy)
  out <- lapply(c(uncorrected=z2raw,corrected=z2corr),function(z2base) matrix(NA_real_,B,4))
  for(bi in seq_len(B)) {
    if(!is.finite(rho_b[bi])||rho_b[bi]<=0) next
    nm <- c(zmean,yz,x1z,x2z)+net_delta[,bi]
    nm <- nm*ratio1[bi]
    om <- ordinary_mean+ordinary_delta[,bi]
    rhs <- c(om[1],nm[2],om[7],om[8])
    for(method in names(out)) {
      z2base <- if(method=="corrected") z2corr else z2raw
      z2b <- (z2base+z2_delta[bi])*ratio2[bi]
      M <- matrix(c(1,nm[1],om[2],om[3],
                    nm[1],z2b,nm[3],nm[4],
                    om[2],nm[3],om[4],om[6],
                    om[3],nm[4],om[6],om[5]),4,4)
      out[[method]][bi,] <- safe_solve(M,rhs)$beta
    }
  }
  base <- lapply(c(uncorrected=z2raw,corrected=z2corr),function(z2base) {
    M <- matrix(c(1,zmean,mean(X[,1]),mean(X[,2]),
                  zmean,z2base,x1z,x2z,
                  mean(X[,1]),x1z,mean(X[,1]^2),mean(X[,1]*X[,2]),
                  mean(X[,2]),x2z,mean(X[,1]*X[,2]),mean(X[,2]^2)),4,4)
    safe_solve(M,c(mean(y),yz,mean(y*X[,1]),mean(y*X[,2])))
  })
  lapply(names(out),function(method) {
    draws <- out[[method]]
    ci <- apply(draws,2,quantile,c(.025,.975),na.rm=TRUE,names=FALSE)
    list(method=method,beta=base[[method]]$beta,lower=ci[1,],upper=ci[2,],
         se=apply(draws,2,sd,na.rm=TRUE),status=base[[method]]$status,
         rcond=base[[method]]$rcond,valid_boot=colSums(is.finite(draws)))
  })
}

result_rows <- function(fit,motif,delta,n,m,rep,method,interval,mean_degree,seeds) {
  data.frame(motif=motif,delta=delta,n=n,m=m,replicate=rep,method=method,
    interval=interval,coefficient=names(truth),estimate=fit$beta,target=unname(truth),
    lower=fit$lower,upper=fit$upper,width=fit$upper-fit$lower,
    covered=unname(truth)>=fit$lower&unname(truth)<=fit$upper,se=fit$se,
    status=fit$status,rcond=fit$rcond,mean_degree=mean_degree,
    data_seed=seeds[1],subset_seed=seeds[2],bootstrap_seed=seeds[3])
}

validate_counts <- function(g,rho_hat,weights) {
  n <- length(g$degree); C2 <- choose(n-1,2); C4 <- choose(n-1,4)
  stopifnot(all(g$triangles>=0),all(g$bowtie_center>=-1e-8),
    abs(sum(g$triangle_squared_inclusion)-5*sum(g$bowtie_center))<1e-6,
    abs(sum(g$two_star_squared_inclusion)-5*sum(choose(g$degree,4)))<1e-6)
  zts <- choose(g$degree,2)/(C2*rho_hat^2)
  ztr <- g$triangles/(C2*rho_hat^3)
  lhs_ts <- colMeans(g$two_star_inclusion/(3*C2*rho_hat^2))
  lhs_tr <- colMeans(g$triangle_inclusion/(C2*rho_hat^3))
  stopifnot(max(abs(lhs_ts-colMeans(weights*zts)))<1e-8,
            max(abs(lhs_tr-colMeans(weights*ztr)))<1e-8)
  qts <- choose(g$degree,4)/(C4*rho_hat^4)
  qtr <- g$bowtie_center/(3*C4*rho_hat^6)
  stopifnot(abs(mean(g$two_star_squared_inclusion/(5*C4*rho_hat^4))-mean(qts))<1e-8,
            abs(mean(g$triangle_squared_inclusion/(15*C4*rho_hat^6))-mean(qtr))<1e-8)
  invisible(TRUE)
}

summarize_results <- function(x) {
  keys <- c("motif","delta","n","m","method","interval","coefficient")
  spl <- split(x,interaction(x[keys],drop=TRUE,lex.order=TRUE))
  ans <- lapply(spl,function(a) {
    ok <- is.finite(a$estimate)&is.finite(a$width)
    data.frame(a[1,keys],repetitions=nrow(a),valid=sum(ok),failure_rate=mean(!ok),
      coverage=if(any(ok))mean(a$covered[ok])else NA_real_,
      coverage_mcse=if(any(ok))sqrt(mean(a$covered[ok])*(1-mean(a$covered[ok]))/sum(ok))else NA_real_,
      mean_width=if(any(ok))mean(a$width[ok])else NA_real_,
      width_mcse=if(sum(ok)>1)sd(a$width[ok])/sqrt(sum(ok))else NA_real_,
      median_width=if(any(ok))median(a$width[ok])else NA_real_,
      mean_estimate=if(any(ok))mean(a$estimate[ok])else NA_real_,
      bias=if(any(ok))mean(a$estimate[ok]-a$target[ok])else NA_real_,
      rmse=if(any(ok))sqrt(mean((a$estimate[ok]-a$target[ok])^2))else NA_real_)
  })
  do.call(rbind,ans)
}

run_experiment <- function(out_dir,n=6000L,repetitions=400L,B=999L,
                           deltas=c(-.55,-.4,-.3,-.25),master_seed=260914L) {
  if(dir.exists(out_dir)) stop("Choose a new output directory")
  dir.create(out_dir,recursive=TRUE)
  m_grid <- main_m_grid[main_m_grid<n]
  config <- list(n=n,repetitions=repetitions,B=B,deltas=deltas,m_grid=m_grid,
    master_seed=master_seed,truth=truth,graphon="w(u,v)=u+v; rho=n^delta",
    response="Y=1+20*Z+3*X1+2*X2+epsilon",
    primary_methods=c("full uncorrected OLS","full bias-corrected OLS","downsampled OLS"),
    diagnostic="bias-corrected quadratic moment within each downsample",
    downsample_interval="exact quantiles of linearized iid Gaussian multiplier distribution",
    full_interval="order-aware first-projection Gaussian multiplier bootstrap",
    session=capture.output(sessionInfo()))
  saveRDS(config,file.path(out_dir,"configuration.rds"))
  rows <- list(); diagnostics <- list(); ri <- 0L
  for(k in seq_along(deltas)) for(rep in seq_len(repetitions)) {
    delta <- deltas[k]; rho <- n^delta
    seeds <- c(master_seed+100000L*k+rep,
               master_seed+1000000L+100000L*k+rep,
               master_seed+2000000L+100000L*k+rep)
    set.seed(seeds[1]); xi <- runif(n)
    X <- sweep(matrix(rnorm(n*2),n,2)%*%chol(matrix(c(1,.6,.6,4),2)),2,c(1,3),"+")
    eps <- rnorm(n); clean <- clean_covariates(xi)
    ys <- lapply(clean,function(z) drop(1+20*z+3*X[,1]+2*X[,2]+eps))
    weights <- cbind(1,ys$two_star,X,1,ys$triangle,X)
    g <- count_subgraphs_cpp(xi,rho,weights)
    rho_hat <- 2*g$edges/(n*(n-1)); if(rho_hat<=0) stop("empty graph")
    if(rep==1L) validate_counts(g,rho_hat,weights)
    C2 <- choose(n-1,2); C4 <- choose(n-1,4)
    observed <- list(two_star=choose(g$degree,2)/(C2*rho_hat^2),
                     triangle=g$triangles/(C2*rho_hat^3))
    corrected <- list(two_star=choose(g$degree,4)/(C4*rho_hat^4),
                      triangle=g$bowtie_center/(3*C4*rho_hat^6))
    set.seed(seeds[2]); permutation <- sample.int(n)
    diagnostics[[length(diagnostics)+1L]] <- data.frame(delta=delta,replicate=rep,
      rho=rho,rho_hat=rho_hat,mean_degree=mean(g$degree),
      zero_triangle_fraction=mean(g$triangles==0),data_seed=seeds[1],subset_seed=seeds[2])
    for(motif in names(clean)) {
      z <- observed[[motif]]; y <- ys[[motif]]; qcorr <- corrected[[motif]]
      inc <- if(motif=="two_star") g$two_star_inclusion[,1:4,drop=FALSE] else g$triangle_inclusion[,5:8,drop=FALSE]
      inc2 <- if(motif=="two_star") g$two_star_squared_inclusion else g$triangle_squared_inclusion
      full <- full_structural_bootstrap(motif,z,y,X,qcorr,inc,inc2,g$degree,rho_hat,B,seeds[3]+match(motif,names(clean)))
      for(f in full) {
        f$lower <- f$lower; f$upper <- f$upper
        ri <- ri+1L; rows[[ri]] <- result_rows(f,motif,delta,n,n,rep,
          paste0("full_",f$method),"structural_multiplier",mean(g$degree),seeds)
      }
      for(m in m_grid) {
        ids <- permutation[seq_len(m)]
        fu <- moment_fit(z,y,X,z^2,ids)
        ri <- ri+1L; rows[[ri]] <- result_rows(fu,motif,delta,n,m,rep,
          "downsampled_ols","iid_gaussian_multiplier",mean(g$degree),seeds)
        fc <- moment_fit(z,y,X,qcorr,ids)
        ri <- ri+1L; rows[[ri]] <- result_rows(fc,motif,delta,n,m,rep,
          "downsampled_corrected_diagnostic","iid_gaussian_multiplier",mean(g$degree),seeds)
      }
    }
    if(rep%%10L==0L||rep==repetitions) message("delta=",delta,": ",rep,"/",repetitions)
    if(rep%%10L==0L||rep==repetitions) {
      current <- do.call(rbind,rows); saveRDS(current,file.path(out_dir,"replicate_results.rds"))
      write.csv(do.call(rbind,diagnostics),file.path(out_dir,"graph_diagnostics.csv"),row.names=FALSE)
    }
  }
  results <- do.call(rbind,rows); summary <- summarize_results(results)
  write.csv(results,file.path(out_dir,"replicate_results.csv"),row.names=FALSE)
  write.csv(summary,file.path(out_dir,"summary_all_coefficients.csv"),row.names=FALSE)
  write.csv(subset(summary,coefficient=="network"),file.path(out_dir,"summary_network.csv"),row.names=FALSE)
  invisible(summary)
}

if(sys.nframe()==0L) {
  a <- commandArgs(trailingOnly=TRUE)
  out <- if(length(a))a[1]else file.path(script_dir,"results","main_n6000_r400")
  n <- if(length(a)>=2)as.integer(a[2])else 6000L
  repetitions <- if(length(a)>=3)as.integer(a[3])else 400L
  B <- if(length(a)>=4)as.integer(a[4])else 999L
  deltas <- if(length(a)>=5)as.numeric(strsplit(a[5],",",fixed=TRUE)[[1]]) else c(-.55,-.4,-.3,-.25)
  run_experiment(out,n,repetitions,B,deltas=deltas)
}
