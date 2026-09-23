args <- commandArgs(trailingOnly=TRUE)
script_arg <- grep("^--file=",commandArgs(),value=TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=","",script_arg)))
root <- if(length(args)) args[1] else file.path(script_dir,"results")
dirs <- file.path(root,c("main_n6000_r400_delta_m055","main_n6000_r400_delta_m040",
                         "main_n6000_r400_delta_m030","main_n6000_r400_delta_m025"))
if(!all(file.exists(file.path(dirs,"replicate_results.rds")))) stop("Not all split runs are available")
source(file.path(script_dir,"coverage_experiment.R"))
x <- do.call(rbind,lapply(dirs,function(d)readRDS(file.path(d,"replicate_results.rds"))))
out <- file.path(root,"main_n6000_r400_combined")
if(!dir.exists(out))dir.create(out,recursive=TRUE)
saveRDS(x,file.path(out,"replicate_results.rds"))
write.csv(x,file.path(out,"replicate_results.csv"),row.names=FALSE)
s <- summarize_results(x)
write.csv(s,file.path(out,"summary_all_coefficients.csv"),row.names=FALSE)
sn <- subset(s,coefficient=="network")
write.csv(sn,file.path(out,"summary_network.csv"),row.names=FALSE)
d <- do.call(rbind,lapply(dirs,function(z)read.csv(file.path(z,"graph_diagnostics.csv"))))
write.csv(d,file.path(out,"graph_diagnostics.csv"),row.names=FALSE)
cfg <- readRDS(file.path(dirs[1],"configuration.rds")); cfg$deltas <- sort(unique(x$delta))
saveRDS(cfg,file.path(out,"configuration.rds"))

cols <- c(downsampled_ols="#2166AC",downsampled_corrected_diagnostic="#D99000",
          full_uncorrected="#B2182B",full_corrected="#1B7837")
pch <- c(downsampled_ols=16,downsampled_corrected_diagnostic=17,
         full_uncorrected=15,full_corrected=18)
labels <- c(downsampled_ols="Downsampled OLS",
            downsampled_corrected_diagnostic="Corrected downsample (diagnostic)",
            full_uncorrected="Full OLS",full_corrected="Full corrected OLS")
for(motif in c("two_star","triangle")) {
  a <- sn[sn$motif==motif,]
  primary <- c("downsampled_ols","full_uncorrected","full_corrected")
  a <- a[a$method %in% primary,]
  png(file.path(out,paste0(motif,"_coverage_width.png")),2400,1300,res=170)
  par(mfrow=c(2,4),mar=c(4.2,4.3,3.1,1),oma=c(5,1,5,1),family="sans")
  for(metric in c("coverage","mean_width")) for(delta in sort(unique(a$delta))) {
    z <- a[a$delta==delta,]
    plot(NA,xlim=range(z$m),log="x",ylim=if(metric=="coverage")c(0,1)else c(0,max(z$mean_width,na.rm=TRUE)*1.08),
         xlab="Regression size m (log scale)",ylab=if(metric=="coverage")"Coverage"else"Average interval width",
         main=sprintf("delta = %.2f | degree %.1f",delta,mean(x$mean_degree[x$delta==delta])))
    if(metric=="coverage")abline(h=.95,lty=2,col="#444444")
    for(method in primary) {
      zz <- z[z$method==method,]; if(!nrow(zz))next
      ord <- order(zz$m); zz <- zz[ord,]
      yy <- zz[[metric]]
      if(grepl("downsampled",method))lines(zz$m,yy,col=cols[method],lwd=1.5,lty=if(method=="downsampled_ols")1 else 3)
      points(zz$m,yy,col=cols[method],pch=pch[method],cex=.8)
      if(metric=="coverage") {
        keep <- is.finite(zz$coverage_mcse)&zz$coverage_mcse>0
        arrows(zz$m[keep],pmax(0,yy[keep]-1.96*zz$coverage_mcse[keep]),
               zz$m[keep],pmin(1,yy[keep]+1.96*zz$coverage_mcse[keep]),
               angle=90,code=3,length=.025,col=adjustcolor(cols[method],.65))
      }
    }
  }
  title_text <- if(motif=="two_star")"Local two-star"else"Local triangle"
  mtext(paste0(title_text,": 95% coverage and interval-width tradeoff"),side=3,outer=TRUE,line=3,cex=1.25,font=2)
  mtext("n = 6000 | 400 independent graphs per sparsity | target network coefficient = 20",side=3,outer=TRUE,line=1.5,cex=.95)
  legend("bottom",inset=-.33,xpd=NA,horiz=TRUE,legend=labels[primary],
         col=cols[primary],pch=pch[primary],lty=c(1,NA,NA),bty="n",cex=.88)
  mtext("Error bars are 95% Monte Carlo intervals for estimated coverage. Full estimators use the structural multiplier; downsampled estimators use the iid Gaussian multiplier.",
        side=1,outer=TRUE,line=3.2,cex=.78)
  dev.off()
}
cat("Wrote combined data, summaries, and figures to",out,"\n")
