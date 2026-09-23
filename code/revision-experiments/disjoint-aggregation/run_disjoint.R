args <- commandArgs(trailingOnly=TRUE)
here <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])))
old_dir <- here
anchor_dir <- normalizePath(file.path(here,"..","estimator-comparison"))
source(file.path(here,"aggregation_core.R"))
source(file.path(here,"disjoint_core.R"))
stopifnot(length(args)>=3)
out <- args[1]; R <- as.integer(args[2]); delta <- as.numeric(args[3])
stopifnot(R>=2L,R<=400L,delta %in% c(-.55,-.4,-.3,-.25))
if(dir.exists(out))stop("Existing results preserved: choose a new output directory")
tag <- sprintf("%03d",as.integer(round(abs(delta)*100)))
anchor_results <- file.path(anchor_dir,"results",paste0("main_n6000_r400_delta_m",tag))
seeds_all <- read.csv(file.path(anchor_results,"graph_diagnostics.csv"))
anchor <- load_anchor(anchor_dir,compile=FALSE)
Rcpp::sourceCpp(file.path(old_dir,"count_covariates.cpp"))
n<-6000L; m_grid<-c(25L,50L,100L); K_grid<-c(1L,5L,10L)
bank<-setNames(lapply(m_grid,function(m)fixed_blocks(n,m,max(K_grid))),m_grid)
dir.create(out,recursive=TRUE)
config<-list(n=n,repetitions=R,delta=delta,m_grid=m_grid,K_grid=K_grid,
  truth=anchor$truth,anchor_directory=anchor_results,
  sampling="Fixed I_k={(k-1)m+1,...,km}; same blocks in every dataset; nested K; no selection RNG",
  feature_estimation="Full graph; sparsity estimated by empirical edge density",
  interval="HC0 independent node Gaussian multiplier; V=sum(V_k)/K^2; analytic normal quantiles, no inner draws",
  methods=c("disjoint_plugin","oracle_diagnostic"),
  caveat="Disjointness does not make graph plug-in estimates exactly independent. Justification uses fixed-K oracle equivalence under downsampling conditions.",
  seed_policy="Replay original saved data seeds; no changes to graph RNG or original source files",
  code_md5=tools::md5sum(c(file.path(here,c("run_disjoint.R","disjoint_core.R")),
    file.path(old_dir,c("aggregation_core.R","count_covariates.cpp")),
    file.path(anchor_dir,"coverage_experiment.R"))),
  rng=RNGkind(),session=capture.output(sessionInfo()))
saveRDS(config,file.path(out,"configuration.rds"))
writeLines(capture.output(str(config)),file.path(out,"configuration.txt"))
saveRDS(bank,file.path(out,"fixed_subset_indices.rds"))
rows<-blocks<-diags<-checks<-list(); ri<-bi<-0L; start<-proc.time()[3]
for(rep in seq_len(R)) {
  seeds<-seeds_all[seeds_all$replicate==rep,]
  stopifnot(nrow(seeds)==1)
  set.seed(seeds$data_seed); xi<-runif(n)
  X<-sweep(matrix(rnorm(n*2),n,2)%*%chol(matrix(c(1,.6,.6,4),2)),2,c(1,3),"+")
  eps<-rnorm(n); clean<-anchor$clean_covariates(xi)
  ys<-lapply(clean,function(z)drop(1+20*z+3*X[,1]+2*X[,2]+eps))
  g<-count_covariates_cpp(xi,n^delta)
  rho_hat<-2*g$edges/(n*(n-1))
  stopifnot(abs(rho_hat-seeds$rho_hat)<1e-14,abs(mean(g$degree)-seeds$mean_degree)<1e-9)
  zs<-list(two_star=choose(g$degree,2)/(choose(n-1,2)*rho_hat^2),
    triangle=g$triangles/(choose(n-1,2)*rho_hat^3))
  for(motif in names(zs)) for(m in m_grid) for(method in config$methods) {
    z<-if(method=="disjoint_plugin")zs[[motif]]else clean[[motif]]
    D<-cbind(1,z,X); sets<-bank[[as.character(m)]]
    a<-disjoint_aggregate(D,ys[[motif]],sets,K_grid)
    # Independent reference fitting code for every K=1 estimate and SE.
    ref<-anchor$moment_fit(z,ys[[motif]],X,q=z^2,ids=sets[[1]])
    dif<-max(abs(c(a$aggregates[[1]]$beta-ref$beta,a$aggregates[[1]]$se-ref$se)))
    stopifnot(dif<1e-8)
    checks[[length(checks)+1L]]<-data.frame(replicate=rep,motif=motif,m=m,method=method,max_error=dif)
    for(k in seq_along(a$blocks)) {
      f<-a$blocks[[k]]; bi<-bi+1L
      blocks[[bi]]<-data.frame(replicate=rep,motif=motif,delta=delta,m=m,method=method,k=k,
        coefficient=names(anchor$truth),estimate=f$beta,variance=colSums(f$contribution^2),
        status=f$status)
    }
    for(f in a$aggregates) {
      ri<-ri+1L
      rows[[ri]]<-data.frame(motif=motif,delta=delta,n=n,m=m,K=f$K,replicate=rep,method=method,
        coefficient=names(anchor$truth),estimate=f$beta,target=unname(anchor$truth),
        se=f$se,lower=f$lower,upper=f$upper,width=f$upper-f$lower,
        covered=f$lower<=anchor$truth & anchor$truth<=f$upper,status=f$status,min_rcond=f$min_rcond,
        unique_nodes=m*f$K,mean_pair_overlap=if(f$K==1)NA_real_ else 0,
        overlap_variance=0,rho_hat=rho_hat,data_seed=seeds$data_seed)
    }
  }
  diags[[rep]]<-data.frame(seeds,elapsed_seconds=unname(proc.time()[3]-start))
  if(rep%%20L==0L || rep==R) {
    saveRDS(do.call(rbind,rows),file.path(out,"replicate_results.rds"))
    saveRDS(do.call(rbind,blocks),file.path(out,"block_results.rds"))
    write.csv(do.call(rbind,diags),file.path(out,"graph_diagnostics.csv"),row.names=FALSE)
    message("delta=",delta,": ",rep,"/",R,"; elapsed ",round(proc.time()[3]-start,1),"s")
  }
}
x<-do.call(rbind,rows)
s<-summarize_aggregation(x); s$expected_pair_overlap<-ifelse(s$K>1,0,NA_real_)
write.csv(s,file.path(out,"summary_all_coefficients.csv"),row.names=FALSE)
write.csv(s[s$coefficient=="network",],file.path(out,"summary_network.csv"),row.names=FALSE)
write.csv(do.call(rbind,checks),file.path(out,"reference_fit_checks.csv"),row.names=FALSE)
old<-readRDS(file.path(anchor_results,"replicate_results.rds"))
old<-old[old$replicate<=R & old$method %in% c("full_uncorrected","full_corrected"),]
saveRDS(old,file.path(out,"paired_full_benchmarks.rds"))
writeLines(paste("Completed",R,"datasets at delta",delta),file.path(out,"COMPLETE.txt"))
