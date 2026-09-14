#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)[1]
repo <- dirname(dirname(normalizePath(sub('^--file=', '', arg))))
out <- file.path(repo,'processed/empirical')
agg <- fread(file.path(out,'option_day_expiry_bins.csv'))
agg[,date:=as.Date(date)]
grid <- seq(-.25,.25,length.out=21)

construct <- function(upper=90, support=TRUE, span=.6) {
  all_dates <- seq(as.Date('2025-01-01'),as.Date('2025-12-31'),by='day')
  iv <- usd <- matrix(NA_real_,365,21)
  details <- vector('list',365)
  for(i in seq_along(all_dates)) {
    dd<-agg[date==all_dates[i] & tau>=7 & tau<=upper]
    taus<-sort(unique(dd$tau));fits<-list()
    for(j in seq_along(taus)) {
      d<-dd[tau==taus[j]]
      if(nrow(d)<8 || (support && (min(d$min_u)>min(grid) || max(d$max_u)<max(grid))))next
      f1<-tryCatch(predict(smooth.spline(d$u_bin,d$iv_dec,spar=span),grid)$y,error=function(e)rep(NA_real_,21))
      f2<-tryCatch(predict(smooth.spline(d$u_bin,log1p(d$usd_notional),spar=span),grid)$y,error=function(e)rep(NA_real_,21))
      if(all(is.finite(c(f1,f2))) && all(f1>0) && all(f2>=0))fits[[as.character(taus[j])]]<-list(iv=f1,usd=f2)
    }
    good<-as.numeric(names(fits));lo<-good[good<30];hi<-good[good>30]
    lower<-upper_used<-NA_real_;mode<-'unavailable';ntr<-0L
    if(length(lo) && length(hi)) {
      lower<-max(lo);upper_used<-min(hi);w<-(30-lower)/(upper_used-lower)
      f1<-fits[[as.character(lower)]];f2<-fits[[as.character(upper_used)]]
      # Total implied variance is linear in tenor; no nearest-tenor fallback.
      iv[i,]<-sqrt(((1-w)*lower*f1$iv^2+w*upper_used*f2$iv^2)/30)
      # This is an interpolated log trading-intensity curve, not traded volume.
      usd[i,]<-(1-w)*f1$usd+w*f2$usd
      mode<-'interpolated'
      ntr<-sum(dd[abs(tau-lower)<1e-8 | abs(tau-upper_used)<1e-8]$n_trades)
      stopifnot(ntr>0)
    }
    details[[i]]<-data.table(date=all_dates[i],mode=mode,lower=lower,upper=upper_used,
                            n_supported_expiries=length(good),selected_trades=ntr)
  }
  list(dates=all_dates,iv=iv,usd=usd,diagnostics=rbindlist(details),grid=grid)
}
main<-construct()
saveRDS(main,file.path(out,'constant30_smiles.rds'))
fwrite(main$diagnostics,file.path(out,'constant30_diagnostics.csv'))
for(ch in c('iv','usd')) {
  curve<-as.data.table(main[[ch]]);setnames(curve,paste0('u',seq_along(grid)))
  fwrite(cbind(data.table(date=main$dates),curve),file.path(out,paste0('constant30_',ch,'_curves.csv')))
}
for(variant in c('upper60','spar05','spar07','extrapolate')) {
  obj<-switch(variant,upper60=construct(upper=60),spar05=construct(span=.5),spar07=construct(span=.7),extrapolate=construct(support=FALSE))
  saveRDS(obj,file.path(out,paste0('constant30_',variant,'.rds')))
  cat(variant,'complete dates',sum(complete.cases(obj$iv)), '\n')
}
spot<-fread(file.path(out,'spot_rv_alternative_sampling.csv'))
spot[,date:=as.Date(date)]
ok<-complete.cases(main$iv)&main$dates %in% spot$date[is.finite(spot$rv_5min)]
r<-rle(ok);ends<-cumsum(r$lengths);starts<-ends-r$lengths+1
runs<-data.table(start=main$dates[starts],end=main$dates[ends],n=r$lengths,complete=r$values)
fwrite(runs,file.path(out,'consecutive_calendar_runs.csv'))
print(runs[complete==TRUE][order(-n)])
cat('Complete fixed-tenor smiles',sum(complete.cases(main$iv)), '\n')
