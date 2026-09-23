script_dir<-dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])))
Rcpp::sourceCpp(file.path(script_dir,"count_covariates.cpp"))
for(seed in 1:4)for(rho in c(.01,.1,.35)) {
  n<-30L;xi<-seq(.01,.99,length.out=n)
  set.seed(seed);g<-count_covariates_cpp(xi,rho);next_cpp<-runif(1)
  set.seed(seed);A<-matrix(0,n,n)
  for(i in 1:(n-1))for(j in (i+1):n) {
    A[i,j]<-A[j,i]<-as.numeric(runif(1)<rho*(xi[i]+xi[j]))
  }
  next_R<-runif(1)
  stopifnot(identical(next_cpp,next_R),identical(g$degree,rowSums(A)),
    identical(g$triangles,diag(A%*%A%*%A)/2),g$edges==sum(A)/2)
}
message("PASS: graph degrees, rooted triangles, edges and RNG continuation match direct adjacency calculations for 12 toy graphs")
