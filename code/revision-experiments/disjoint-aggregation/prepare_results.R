here<-dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])))
source(file.path(here,"disjoint_core.R"))
tags<-c("055","040","030","025"); deltas<-c(-.55,-.4,-.3,-.25)
dirs<-file.path(here,"results",paste0("formal_n6000_r400_delta_m",tags))
stopifnot(all(file.exists(file.path(dirs,"COMPLETE.txt"))))
out<-file.path(here,"results","formal_n6000_r400_combined")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
key<-function(x,cols) do.call(paste,c(x[cols],sep="|"))
near<-function(a,b,tol=1e-8) stopifnot(length(a)==length(b),all(is.finite(a)),all(is.finite(b)),max(abs(a-b))<tol)
keys<-c("motif","delta","n","m","K","coefficient","method")
summaries<-records<-audits<-covdiags<-paired<-benchmarks<-list()
for(j in seq_along(dirs)) {
  dir<-dirs[j]; cfg<-readRDS(file.path(dir,"configuration.rds"))
  x<-readRDS(file.path(dir,"replicate_results.rds")); b<-readRDS(file.path(dir,"block_results.rds"))
  s<-read.csv(file.path(dir,"summary_all_coefficients.csv"))
  stopifnot(nrow(x)==400*2*3*3*2*4,nrow(b)==400*2*3*2*10*4,
    nrow(s)==2*3*3*2*4,all(s$valid==400),all(s$failed==0),all(x$status=="ok"),all(b$status=="ok"),
    !anyDuplicated(key(x,c(keys,"replicate"))),all(x$unique_nodes==x$m*x$K),
    all(x$overlap_variance==0),all(x$mean_pair_overlap[x$K>1]==0),
    identical(unname(tools::md5sum(names(cfg$code_md5))),unname(cfg$code_md5)))
  bank<-readRDS(file.path(dir,"fixed_subset_indices.rds"))
  for(m in cfg$m_grid) stopifnot(identical(unlist(bank[[as.character(m)]]),seq_len(10L*m)))
  # Reconstruct the reported aggregates directly from the stored separate fits.
  join<-c("replicate","motif","delta","m","method","coefficient")
  for(K in cfg$K_grid) {
    bb<-b[b$k<=K,]; group<-key(bb,join); xx<-x[x$K==K,]; q<-key(xx,join)
    est<-tapply(bb$estimate,group,mean)
    v<-tapply(bb$variance,group,sum)/K^2
    near(xx$estimate,as.numeric(est[q]));near(xx$se^2,as.numeric(v[q]))
  }
  g<-key(x,keys); q<-key(s,keys); stopifnot(all(table(g)==400))
  hit<-as.numeric(x$lower<=x$target & x$target<=x$upper)
  near(x$covered,hit);near(x$lower,x$estimate-qnorm(.975)*x$se)
  near(x$upper,x$estimate+qnorm(.975)*x$se)
  near(s$coverage,as.numeric(tapply(hit,g,mean)[q]))
  near(s$mean_width,as.numeric(tapply(x$upper-x$lower,g,mean)[q]))
  near(s$bias,as.numeric(tapply(x$estimate-x$target,g,mean)[q]))
  near(s$empirical_sd,as.numeric(tapply(x$estimate,g,sd)[q]))
  near(s$rmse,sqrt(as.numeric(tapply((x$estimate-x$target)^2,g,mean)[q])))
  near(s$rms_se,sqrt(as.numeric(tapply(x$se^2,g,mean)[q])))
  z<-qnorm(.975); denom<-1+z^2/s$valid
  ctr<-(s$coverage+z^2/(2*s$valid))/denom
  half<-z*sqrt(s$coverage*(1-s$coverage)/s$valid+z^2/(4*s$valid^2))/denom
  s$coverage_mc_lower<-pmax(0,ctr-half);s$coverage_mc_upper<-pmin(1,ctr+half)
  s$se_to_empirical_sd<-s$rms_se/s$empirical_sd
  summaries[[j]]<-s;records[[j]]<-x[x$coefficient=="network",]
  # Empirical variance decomposition of the network coefficient across datasets.
  for(motif in c("two_star","triangle")) for(m in cfg$m_grid) for(method in cfg$methods) {
    bb<-b[b$motif==motif & b$m==m & b$method==method & b$coefficient=="network",]
    # Construct robustly by explicit filtering, avoiding any dependence on row layout.
    B<-sapply(1:10,function(k) {a<-bb[bb$k==k,];a$estimate[match(1:400,a$replicate)]})
    for(K in cfg$K_grid) {
      C<-cov(B[,seq_len(K),drop=FALSE]); total<-var(rowMeans(B[,seq_len(K),drop=FALSE]))
      diagonal<-sum(diag(C))/K^2; cross<-sum(C[row(C)!=col(C)])/K^2
      near(total,diagonal+cross)
      aa<-x[x$motif==motif & x$m==m & x$method==method & x$coefficient=="network" & x$K==K,]
      near(aa$estimate[match(1:400,aa$replicate)],rowMeans(B[,seq_len(K),drop=FALSE]))
      covdiags[[length(covdiags)+1L]]<-data.frame(motif=motif,delta=cfg$delta,m=m,method=method,K=K,
        empirical_variance=total,diagonal_variance=diagonal,cross_covariance_contribution=cross,
        rms_se=sqrt(mean(aa$se^2)),variance_ratio_to_K1=total/var(B[,1]),
        independent_reference_ratio=1/K)
      if(K>1) {
        base<-x[x$motif==motif & x$m==m & x$method==method & x$coefficient=="network" & x$K==1,]
        base<-base[match(aa$replicate,base$replicate),]
        de<-aa$estimate-base$estimate; dc<-as.numeric(aa$covered)-as.numeric(base$covered)
        paired[[length(paired)+1L]]<-data.frame(motif=motif,delta=cfg$delta,m=m,method=method,K=K,
          estimate_change=mean(de),estimate_change_mcse=sd(de)/20,
          coverage_change=mean(dc),coverage_change_mcse=sd(dc)/20)
      }
    }
  }
  checks<-read.csv(file.path(dir,"reference_fit_checks.csv"));stopifnot(max(checks$max_error)<1e-8)
  audits[[j]]<-data.frame(delta=cfg$delta,replications=400,coefficient_records=nrow(x),cells=nrow(s),
    failures=0,disjointness_verified=TRUE,aggregate_variance_reconstructed=TRUE,
    summaries_recomputed=TRUE,code_hashes_match=TRUE,max_reference_fit_error=max(checks$max_error))
  benchmarks[[j]]<-readRDS(file.path(dir,"paired_full_benchmarks.rds"))
}
s<-do.call(rbind,summaries); x<-do.call(rbind,records)
write.csv(s,file.path(out,"summary_all_coefficients.csv"),row.names=FALSE)
net<-s[s$coefficient=="network",]
write.csv(net,file.path(out,"summary_network.csv"),row.names=FALSE)
saveRDS(x,file.path(out,"replicate_network.rds"))
write.csv(do.call(rbind,audits),file.path(out,"validation_checks.csv"),row.names=FALSE)
write.csv(do.call(rbind,covdiags),file.path(out,"variance_diagnostics.csv"),row.names=FALSE)
write.csv(do.call(rbind,paired),file.path(out,"paired_changes_vs_K1.csv"),row.names=FALSE)
bs<-do.call(rbind,benchmarks); bkeys<-c("motif","delta","method","coefficient")
bs<-do.call(rbind,lapply(split(bs,key(bs,bkeys)),function(b) data.frame(b[1,bkeys],
  coverage=mean(b$covered),mean_width=mean(b$width),bias=mean(b$estimate-b$target),
  empirical_sd=sd(b$estimate),valid=nrow(b))))
write.csv(bs,file.path(out,"full_benchmark_summary.csv"),row.names=FALSE)

primary<-net[net$method=="disjoint_plugin",];m_grid<-c(25,50,100);K_grid<-c(1,5,10)
cell<-function(p,w)sprintf("%.3f (%.2f)",floor(p*1000+.5+1e-9)/1000,w)
tex<-c("% Requires booktabs. Fixed disjoint subsets; checked Monte Carlo results.")
for(motif in c("two_star","triangle")) {
  title<-if(motif=="two_star")"local two-star"else"local triangle"
  t<-c("\\begin{table}[tbp]","\\centering","\\small",
    paste0("\\caption{Coverage and interval width for disjoint-subset aggregation with a ",title," covariate.}"),
    paste0("\\label{tab:disjoint-aggregation-",gsub("_","-",motif),"}"),
    "\\begin{tabular}{@{}rrrrrr@{}}","\\toprule",
    "$m$ & $K$ & $\\delta=-0.55$ & $-0.40$ & $-0.30$ & $-0.25$ \\\\","\\midrule")
  for(m in m_grid) {
    for(K in K_grid) {
      vals<-vapply(deltas,function(d) {
        a<-primary[primary$motif==motif & primary$m==m & primary$K==K & primary$delta==d,]
        stopifnot(nrow(a)==1);cell(a$coverage,a$mean_width)
      },character(1))
      t<-c(t,paste0(paste(c(if(K==1)m else "",K,vals),collapse=" & ")," \\\\"))
    }
    if(m!=100)t<-c(t,"\\addlinespace")
  }
  t<-c(t,"\\bottomrule","\\end{tabular}","\\par\\smallskip",
    "{\\footnotesize Entries are coverage of nominal 95\\% intervals (mean width); $n=6000$, 400 replications per regime.}","\\end{table}")
  writeLines(t,file.path(out,paste0("disjoint_",motif,"_table.tex")));tex<-c(tex,t,"")
  a<-primary[primary$motif==motif,]
  draw<-function() {
    par(mfrow=c(3,4),mar=c(3,3.5,2.1,.6),oma=c(3.6,.7,4,.5),
      family="serif",mgp=c(2.2,.65,0),tcl=-.25,cex=.85)
    for(metric in c("coverage","mean_width","bias"))for(d in deltas) {
      b<-a[a$delta==d,]
      yl<-switch(metric,coverage=c(0,1),mean_width=c(0,max(a$mean_width+1.96*a$width_mcse)),
        bias=range(c(0,a$bias-1.96*a$bias_mcse,a$bias+1.96*a$bias_mcse)))
      lab<-switch(metric,coverage="Coverage",mean_width="Mean CI width",bias="Bias")
      plot(NA,xlim=c(.4,10.6),ylim=yl,xlab="K",ylab=lab,xaxt="n",bty="l",main=bquote(delta==.(d)))
      axis(1,at=K_grid);abline(h=axTicks(2),col="grey92",lwd=.6)
      if(metric=="coverage")abline(h=.95,lty=3,col="grey45")
      if(metric=="bias")abline(h=0,lty=3,col="grey45")
      for(j in seq_along(m_grid)) {
        c<-b[b$m==m_grid[j],];c<-c[order(c$K),]
        lo<-switch(metric,coverage=c$coverage_mc_lower,mean_width=c$mean_width-1.96*c$width_mcse,
          bias=c$bias-1.96*c$bias_mcse)
        hi<-switch(metric,coverage=c$coverage_mc_upper,mean_width=c$mean_width+1.96*c$width_mcse,
          bias=c$bias+1.96*c$bias_mcse)
        color<-c("#222222","#626262","#929292")[j]
        segments(c$K,lo,c$K,hi,col=color,lwd=.7)
        segments(c$K-.10,lo,c$K+.10,lo,col=color,lwd=.7)
        segments(c$K-.10,hi,c$K+.10,hi,col=color,lwd=.7)
        lines(c$K,c[[metric]],type="b",pch=c(16,17,0)[j],lty=j,col=color,lwd=1.1,cex=.65)
      }
    }
    par(fig=c(0,1,0,1),new=TRUE,mar=rep(0,4),oma=rep(0,4));plot.new()
    plot.window(xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i")
    text(.5,.976,paste0(if(motif=="two_star")"Local two-star"else"Local triangle",": disjoint-subset aggregation"),font=2,cex=1.1)
    legend(x=.5,y=.94,xjust=.5,yjust=.5,legend=paste("m =",m_grid),
      col=c("#222222","#626262","#929292"),pch=c(16,17,0),lty=1:3,horiz=TRUE,bty="n",cex=.85)
    text(.5,.016,"n = 6000; 400 replications per regime. Bars: 95% Monte Carlo uncertainty.",cex=.75)
  }
  png(file.path(out,paste0("disjoint_aggregation_",motif,".png")),width=3000,height=2400,res=300)
  draw();dev.off()
}
writeLines(tex,file.path(out,"disjoint_aggregation_tables.tex"))
writeLines(c("PASS: 230400 coefficient records, 576 summary cells, each with 400 valid replications.",
  "All fixed blocks checked; all estimates and multiplier variances reconstructed from block records.",
  "All reported coverage, width, bias, SD, RMSE and RMS SE independently recomputed.",
  "Code hashes checked; no fit failures. Figures require visual inspection."),file.path(out,"validation_summary.txt"))
message("Prepared and numerically verified: ",out)
