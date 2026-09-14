#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)[1]
repo<-dirname(dirname(normalizePath(sub('^--file=','',arg))))
out<-file.path(repo,'processed/empirical')
source(file.path(repo,'R/FGC_cv_lambda20_fixed.R'))
cv<-c(quad=get_cv_lambda20(.05,'quad'),range=get_cv_lambda20(.05,'range'))
lambda<- .25+.75*(1:19)/20
spot<-fread(file.path(out,'spot_rv_alternative_sampling.csv'));spot[,date:=as.Date(date)]
smiles<-readRDS(file.path(out,'constant30_smiles.rds'))
calendar<-seq(as.Date('2025-01-01'),as.Date('2025-12-31'),by='day')
zscore<-function(A) {
  A<-as.matrix(A);s<-apply(A,2,sd)
  if(any(!is.finite(s)|s<=0))stop('Degenerate coordinate')
  sweep(sweep(A,2,colMeans(A),'-'),2,s,'/')
}
zdesign<-function(Z)cbind(1,zscore(Z[,-1,drop=FALSE]))
resid_ridge<-function(Z,A,mult=1) {
  q<-nrow(Z);pen<-diag(ncol(Z));pen[1,1]<-0
  A-Z%*%solve(crossprod(Z)+mult*pen,crossprod(Z,A))
}
make_rows<-function(obj,channel,direction,rv='spot_vol_5min') {
  y<-spot[[rv]][match(calendar,spot$date)]
  X<-obj[[tolower(channel)]][match(calendar,obj$dates),,drop=FALSE]
  target<-23:365
  if(direction=='Options to RV') {
    Y<-matrix(y[target],ncol=1);XL<-X[target-1,,drop=FALSE]
    Z<-cbind(1,y[target-1],vapply(target,function(t)mean(y[(t-5):(t-1)]),0.),vapply(target,function(t)mean(y[(t-22):(t-1)]),0.))
  } else {
    Y<-X[target,,drop=FALSE];XL<-matrix(y[target-1],ncol=1)
    Z<-cbind(1,X[target-1,,drop=FALSE])
  }
  valid<-complete.cases(Y,XL,Z)
  list(Y=Y,X=XL,Z=Z,valid=valid,dates=calendar[target],T=length(target),Dx=ncol(XL),Dy=ncol(Y),channel=channel,direction=direction)
}
subset_rows<-function(rows,ix) {
  ans<-rows
  for(k in c('Y','X','Z'))ans[[k]]<-rows[[k]][ix,,drop=FALSE]
  ans$valid<-rows$valid[ix];ans$dates<-rows$dates[ix];ans$T<-length(ix);ans
}
fit_moments<-function(rows,ix,mult=1) {
  ix<-ix[rows$valid[ix]]
  if(length(ix)<ncol(rows$Z)+3)stop('Insufficient usable pairs for fixed prefix')
  Y<-zscore(rows$Y[ix,,drop=FALSE]);X<-zscore(rows$X[ix,,drop=FALSE]);Z<-zdesign(rows$Z[ix,,drop=FALSE])
  U<-resid_ridge(Z,Y,mult);V<-resid_ridge(Z,X,mult)
  C<-crossprod(U,V)/length(ix)
  raw_d<-sum(C*C);q<-length(ix)
  diagonal<-mean(rowSums(U^2)*rowSums(V^2))
  # Remove self-products from the squared cross-moment. This correction is
  # O_p(1/q) in fixed dimension, so it preserves the first-order SN expansion.
  # It is not claimed to be exactly unbiased with serially dependent residuals.
  corrected<-(q*raw_d-diagonal)/(q-1)
  list(d=corrected,d_raw=raw_d,C=C,U=U,V=V,n=length(ix),Z=Z)
}
# Calendar fractions, not compressed observation indices. The ratio moment
# estimator uses available next-day pairs within each fixed calendar prefix.
sn<-function(rows,trim=.25,mult=1) {
  lam<-trim+(1-trim)*(1:19)/20;full<-fit_moments(rows,1:rows$T,mult)
  ks<-floor(rows$T*lam)
  fitted<-lapply(ks,function(k)fit_moments(rows,seq_len(k),mult))
  # k/T versus lambda differs by O(1/T), and is fixed independently of I_t.
  path<-(ks/rows$T)^2*vapply(fitted,`[[`,0.,'d')
  G<-path-lam^2*full$d
  c(d=full$d,V=sqrt(mean(G*G)),H=diff(range(G)),N=full$n,
    min_prefix_N=min(vapply(fitted,`[[`,0L,'n')),T=rows$T)
}
result_rows<-function(rows,trim=.25,mult=1) {
  a<-sn(rows,trim,mult)
  r<-data.table(channel=rows$channel,direction=rows$direction,s0=c(.01,.025,.05),
                N=a['N'],T=a['T'],min_prefix_N=a['min_prefix_N'],d_hat=a['d'],s_hat=sqrt(max(a['d'],0)/(rows$Dx*rows$Dy)),V=a['V'],H=a['H'])
  r[,Delta:=rows$Dx*rows$Dy*s0^2]
  r[,`:=`(stat_quad=(d_hat-Delta)/V,stat_range=(d_hat-Delta)/H)]
  r[,`:=`(reject_quad=stat_quad>cv['quad'],reject_range=stat_range>cv['range'])]
  r
}
full<-list();designs<-list();calibration<-list()
for(ch in c('IV','USD'))for(direction in c('Options to RV','RV to Options')) {
  rows<-make_rows(smiles,ch,direction);key<-paste(ch,direction)
  designs[[key]]<-rows;full[[key]]<-result_rows(rows)
  f<-fit_moments(rows,seq_len(rows$T));S<-crossprod(f$V)/f$n
  ev<-eigen(S,symmetric=TRUE);lmax<-max(ev$values)
  use<-ev$values>max(ev$values)*1e-10
  sinv<-ev$vectors[,use,drop=FALSE]%*%diag(1/ev$values[use],sum(use))%*%t(ev$vectors[,use,drop=FALSE])
  gain<-sum((f$C%*%sinv)*f$C);base<-sum(f$U*f$U)/f$n
  calibration[[key]]<-data.table(channel=ch,direction=direction,s0=c(.01,.025,.05),
    guaranteed_population_gain_pct=100*c(.01,.025,.05)^2,
    residual_predictor_lambda_max=lmax,restricted_MSE=base,
    plugin_spectral_floor_pct=100*rows$Dx*rows$Dy*c(.01,.025,.05)^2/(lmax*base),
    sample_projection_gain_pct=100*gain/base)
}
full<-rbindlist(full);fwrite(full,file.path(out,'calendar_full_sample.csv'))
fwrite(rbindlist(calibration),file.path(out,'threshold_gain_calibration.csv'))
saveRDS(designs,file.path(out,'calendar_designs.rds'))

robust<-list()
for(rv in c('1min','5min','10min','15min','5min_subsampled'))for(ch in c('IV','USD')) {
  rows<-make_rows(smiles,ch,'Options to RV',paste0('spot_vol_',rv))
  x<-result_rows(rows)[s0==.05];x[,specification:=rv];robust[[paste(rv,ch)]]<-x
}
fwrite(rbindlist(robust),file.path(out,'calendar_rv_robustness.csv'))
tenor<-list()
for(v in c('upper60','spar05','spar07','extrapolate'))for(ch in c('IV','USD')) {
  obj<-readRDS(file.path(out,paste0('constant30_',v,'.rds')))
  rows<-make_rows(obj,ch,'Options to RV')
  x<-tryCatch(result_rows(rows)[s0==.05],error=function(e)NULL)
  if(!is.null(x)){x[,specification:=v];tenor[[paste(v,ch)]]<-x}
}
fwrite(rbindlist(tenor),file.path(out,'calendar_smile_robustness.csv'))
tenor_controls<-list()
for(ch in c('IV','USD')) {
  rows<-make_rows(smiles,ch,'Options to RV')
  info<-smiles$diagnostics[match(rows$dates-1,smiles$diagnostics$date)]
  rows$Z<-cbind(rows$Z,info$lower,info$upper,log1p(info$selected_trades))
  rows$valid<-rows$valid & complete.cases(rows$Z)
  x<-result_rows(rows)[s0==.05];x[,specification:='tenor_and_trade_count_controls'];tenor_controls[[ch]]<-x
}
fwrite(rbindlist(c(tenor,tenor_controls)),file.path(out,'calendar_smile_robustness.csv'))

# Rolling windows remain in calendar time. Omit an entire window if any of
# its fixed prefixes has too few usable pairs; never drop individual prefixes.
roll<-list();excluded<-list();jj<-kk<-1L
for(key in names(designs)) {
  r<-designs[[key]]
  for(end in seq(180L,r$T,by=5L)) {
    w<-subset_rows(r,(end-179L):end)
    x<-tryCatch(result_rows(w),error=function(e)e)
    if(inherits(x,'error')) {
      excluded[[kk]]<-data.table(channel=r$channel,direction=r$direction,start=w$dates[1],end=w$dates[180],reason=conditionMessage(x));kk<-kk+1L
    } else {
      x[,`:=`(start=w$dates[1],end=w$dates[180])];roll[[jj]]<-x;jj<-jj+1L
    }
  }
}
roll<-rbindlist(roll)
set.seed(20260917)
lam<-lambda;K<-outer(lam,lam,pmin)-outer(lam,lam,'*')
U<-matrix(rnorm(200000*19),200000,19)%*%chol(K)
U<-sweep(U,2,lam,'*');vq<-sqrt(rowMeans(U*U));vr<-apply(U,1,function(x)diff(range(x)))
piv_p<-function(x,v)mean(pnorm(-x*v))
sens_files<-c('calendar_rv_robustness.csv','calendar_smile_robustness.csv')
sens<-rbindlist(lapply(sens_files,function(nm){a<-fread(file.path(out,nm));a[,source_file:=nm];a}))
sens[,p_quad:=vapply(stat_quad,piv_p,0.,v=vq)]
sens[,p_range:=vapply(stat_range,piv_p,0.,v=vr)]
ps<-p.adjust(c(sens$p_quad,sens$p_range),'holm')
sens[,`:=`(p_holm_quad=ps[seq_len(.N)],p_holm_range=ps[nrow(sens)+seq_len(.N)])]
for(nm in sens_files)fwrite(sens[source_file==nm,!c('source_file')],file.path(out,nm))
full[,p_quad:=vapply(stat_quad,piv_p,0.,v=vq)]
full[,p_range:=vapply(stat_range,piv_p,0.,v=vr)]
# The manuscript reports pointwise p-values at three positive thresholds.
# The 0.01 and 0.025 thresholds were added after the original analysis.
# Retain per-threshold Holm comparisons as supplementary diagnostics.
for(s in unique(full$s0)) {
  ix<-which(full$s0==s);a<-p.adjust(c(full$p_quad[ix],full$p_range[ix]),'holm')
  full[ix,`:=`(p_holm_quad=a[seq_along(ix)],p_holm_range=a[length(ix)+seq_along(ix)])]
}
fwrite(full,file.path(out,'calendar_full_sample.csv'))
roll[,p_quad:=vapply(stat_quad,piv_p,0.,v=vq)]
roll[,p_range:=vapply(stat_range,piv_p,0.,v=vr)]
# A single family includes both methods, channels, directions and all windows
# at each prespecified threshold; use the smallest threshold for summaries.
for(s in unique(roll$s0)) {
  ix<-which(roll$s0==s);a<-p.adjust(c(roll$p_quad[ix],roll$p_range[ix]),'holm')
  roll[ix,`:=`(p_holm_quad=a[seq_along(ix)],p_holm_range=a[length(ix)+seq_along(ix)])]
}
fwrite(roll,file.path(out,'calendar_rolling.csv'))
fwrite(rbindlist(excluded),file.path(out,'calendar_rolling_excluded.csv'))

# Forecasts use only earlier calendar dates. The principal comparison predicts
# volatility directly with the same HAR controls and all 21 smile coordinates.
# QLIKE is a separate variance forecast exercise using log-variance HAR and a
# training-residual smearing adjustment. No full-sample PCA or scaling is used.
predict_fit<-function(A,y,a) {
  A<-as.matrix(A);a<-as.numeric(a);mu<-colMeans(A);s<-apply(A,2,sd)
  if(any(s<=0))stop('Degenerate training regressor')
  Z<-cbind(1,sweep(sweep(A,2,mu,'-'),2,s,'/'))
  z<-c(1,(a-mu)/s);pen<-diag(ncol(Z));pen[1,1]<-0
  beta<-solve(crossprod(Z)+pen,crossprod(Z,y))
  list(pred=sum(z*beta),resid=y-as.numeric(Z%*%beta))
}
forecasts<-list();ii<-1L
for(ch in c('IV','USD')) {
  r<-designs[[paste(ch,'Options to RV')]]
  # Log-variance controls are constructed on the complete spot calendar.
  logrv<-log(spot$rv_5min[match(calendar,spot$date)])
  ti<-23:365
  LZ<-cbind(logrv[ti-1],vapply(ti,function(t)mean(logrv[(t-5):(t-1)]),0.),vapply(ti,function(t)mean(logrv[(t-22):(t-1)]),0.))
  for(t in seq_len(r$T)) {
    train<-which(r$valid & seq_len(r$T)<t & seq_len(r$T)>=t-120)
    if(!r$valid[t] || length(train)<60 || t<=120)next
    for(model in c('HAR','HAR+smile')) {
      A<-r$Z[,-1,drop=FALSE];LA<-LZ
      if(model=='HAR+smile'){A<-cbind(A,r$X);LA<-cbind(LA,r$X)}
      f<-predict_fit(A[train,,drop=FALSE],r$Y[train,1],A[t,])
      fl<-predict_fit(LA[train,,drop=FALSE],logrv[ti[train]],LA[t,])
      forecasts[[ii]]<-data.table(channel=ch,model=model,date=r$dates[t],n_train=length(train),
                                  actual_vol=r$Y[t,1],pred_vol=f$pred,actual_variance=r$Y[t,1]^2,
                                  pred_variance=exp(fl$pred)*mean(exp(fl$resid)))
      ii<-ii+1L
    }
  }
}
forecasts<-rbindlist(forecasts);fwrite(forecasts,file.path(out,'calendar_forecasts.csv'))
# Calendar-lag HAC: absent forecast dates contribute no pair at that lag.
hac<-function(d,dates,L) {
  z<-d-mean(d);n<-length(d);v<-sum(z*z)/n
  for(h in seq_len(L)) {
    j<-match(dates-h,dates);ok<-!is.na(j)
    v<-v+2*(1-h/(L+1))*sum(z[ok]*z[j[ok]])/n
  }
  if(!is.finite(v)||v<=0)return(NA_real_)
  2*pnorm(-abs(mean(d)/sqrt(v/n)))
}
loss_results<-list()
for(ch in c('IV','USD')) {
  a<-forecasts[channel==ch & model=='HAR'];b<-forecasts[channel==ch & model=='HAR+smile'];stopifnot(identical(a$date,b$date))
  se0<-(a$actual_vol-a$pred_vol)^2;se1<-(b$actual_vol-b$pred_vol)^2
  qlike<-function(y,h){r<-y/h;r-log(r)-1}
  q0<-qlike(a$actual_variance,a$pred_variance);q1<-qlike(b$actual_variance,b$pred_variance)
  L<-max(1L,floor(nrow(a)^(1/3)))
  x<-data.table(channel=ch,N=nrow(a),rmse_base=sqrt(mean(se0)),rmse_smile=sqrt(mean(se1)),MSE_gain_pct=100*(1-mean(se1)/mean(se0)),
                qlike_base=mean(q0),qlike_smile=mean(q1),QLIKE_gain_pct=100*(1-mean(q1)/mean(q0)),bandwidth=L,
                p_squared_error=hac(se0-se1,a$date,L),p_qlike=hac(q0-q1,a$date,L),
                negative_vol_forecasts=sum(c(a$pred_vol,b$pred_vol)<0),min_training_pairs=min(c(a$n_train,b$n_train)))
  for(l in c(1,3,5,10)){x[[paste0('p_SE_L',l)]]<-hac(se0-se1,a$date,l);x[[paste0('p_QLIKE_L',l)]]<-hac(q0-q1,a$date,l)}
  loss_results[[ch]]<-x
}
fwrite(rbindlist(loss_results),file.path(out,'calendar_forecast_summary.csv'))
print(full);print(rbindlist(loss_results));cat('Calendar analysis complete\n')
