args <- commandArgs(trailingOnly=TRUE)
p <- if(length(args))args[1]else "results/main_n6000_r400_combined"
x <- readRDS(file.path(p,"replicate_results.rds"))
s <- read.csv(file.path(p,"summary_all_coefficients.csv"),check.names=FALSE)
cfg <- readRDS(file.path(p,"configuration.rds"))
expected_fits <- length(cfg$deltas)*cfg$repetitions*2*(2+2*length(cfg$m_grid))
stopifnot(nrow(x)==4*expected_fits,
          all(x$target==unname(c(intercept=1,network=20,x1=3,x2=2)[x$coefficient])),
          all(x$width>=0|is.na(x$width)),
          all(x$covered==(x$target>=x$lower&x$target<=x$upper)|is.na(x$covered)))
cell <- table(x$motif,x$delta,x$m,x$method,x$coefficient)
stopifnot(all(cell[cell>0]==cfg$repetitions))
keys <- c("motif","delta","n","m","method","interval","coefficient")
for(i in seq_len(nrow(s))) {
  a <- x
  for(k in keys)a <- a[a[[k]]==s[[k]][i],,drop=FALSE]
  ok <- is.finite(a$estimate)&is.finite(a$width)
  stopifnot(nrow(a)==cfg$repetitions,sum(ok)==s$valid[i])
  if(any(ok)) stopifnot(abs(mean(a$covered[ok])-s$coverage[i])<1e-12,
                         abs(mean(a$width[ok])-s$mean_width[i])<1e-10,
                         abs(mean(a$estimate[ok])-s$mean_estimate[i])<1e-10,
                         abs(sqrt(mean((a$estimate[ok]-a$target[ok])^2))-s$rmse[i])<1e-10)
}
net <- x[x$coefficient=="network",]
full <- aggregate(cbind(empirical_sd=estimate,mean_bootstrap_se=se)~motif+delta+method,
                  net[grepl("^full",net$method),],
                  function(v)c(sd=sd(v),mean=mean(v)))
write.csv(full,file.path(p,"full_interval_scale_check.csv"),row.names=FALSE)
cat("PASS:",nrow(x),"rows;",cfg$repetitions,"replications in every populated cell;",
    "targets, interval logic, coverage, widths, means, and RMSE all independently verified.\n")
