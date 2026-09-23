here<-dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])))
source(file.path(here,"aggregation_core.R"))
source(file.path(here,"disjoint_core.R"))
near<-function(a,b,tol=1e-10)stopifnot(max(abs(a-b))<tol)
set.seed(20260919); n<-1200L
D<-cbind(1,rnorm(n),rnorm(n),rnorm(n)); y<-drop(D%*%c(1,20,3,2)+rnorm(n))
for(m in c(25L,50L,100L)) {
  sets<-fixed_blocks(n,m,10L)
  stopifnot(identical(unlist(sets),seq_len(10L*m)))
  a<-disjoint_aggregate(D,y,sets)
  shared<-aggregate_subsets(D,y,sets,c(1L,5L,10L))
  for(f in a$aggregates) {
    # Independent lm() fit and direct sandwich calculation, not ols_influence().
    ref<-lapply(sets[seq_len(f$K)],function(ids) {
      d<-D[ids,]; z<-lm.fit(d,y[ids]); inv<-solve(crossprod(d))
      list(beta=z$coefficients,V=inv%*%crossprod(d*z$residuals)%*%inv)
    })
    near(f$beta,Reduce(`+`,lapply(ref,`[[`,"beta"))/f$K)
    near(f$covariance,Reduce(`+`,lapply(ref,`[[`,"V"))/f$K^2)
    near(f$covariance,shared[[as.character(f$K)]]$covariance)
  }
}
# Direct Gaussian draws validate the full covariance, including off-diagonals.
f<-a$aggregates[[3]]; set.seed(3192026)
draws<-matrix(0,40000,4)
for(k in 1:10) draws<-draws+matrix(rnorm(40000*100),40000,100)%*%a$blocks[[k]]$contribution/10
V<-cov(draws); scale<-sqrt(outer(diag(f$covariance),diag(f$covariance)))
stopifnot(max(abs(V-f$covariance)/scale)<.04)
# Invalid allocation must fail; singular fit must not be dropped silently.
stopifnot(inherits(try(disjoint_aggregate(D,y,list(1:25,1:25),2),silent=TRUE),"try-error"))
bad<-disjoint_aggregate(cbind(1,0,D[,3:4]),y,fixed_blocks(n,25),c(1,5,10))
stopifnot(all(vapply(bad$aggregates,function(f)f$status=="subset_fit_failure",logical(1))))
message("PASS: disjoint blocks, lm/sandwich agreement, Gaussian covariance, K=1 equivalence, invalid-allocation and singular-fit handling")
