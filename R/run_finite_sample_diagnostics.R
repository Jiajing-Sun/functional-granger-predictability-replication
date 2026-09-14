#!/usr/bin/env Rscript
source('FGC_simulation_functions_lambda20_plugin.R')
library(data.table)
out<-'../processed/diagnostics';dir.create(out,recursive=TRUE,showWarnings=FALSE)
cvtab<-fread('../processed/simulation/critical_values.csv')
R<-1000L
pop<-function(D,scaled,kind) {
  a<-duX_white_noise_closed_form(.3,.4,.5,1/(1:D),1/(1:D))$components
  c<-a$partial_cov
  if(scaled)c<-c/sqrt(a$var_x*a$var_y)
  if(kind %in% c('scalar','reverse'))sum(c[1]^2) else sum(c*c)
}
z<-function(M)sweep(sweep(M,2,colMeans(M),'-'),2,apply(M,2,sd),'/')
path<-function(dat,scaled,kind,trim=.25,ridge=1) {
  n<-nrow(dat$Y);lam<-c(trim+(1-trim)*(1:19)/20,1)
  vals<-sapply(lam,function(l) {
    k<-floor(n*l);Y<-dat$Y[2:k,,drop=FALSE];X<-dat$X[1:(k-1),,drop=FALSE];Z<-dat$Y[1:(k-1),,drop=FALSE]
    if(kind=='scalar'){Y<-Y[,1,drop=FALSE];Z<-Z[,1,drop=FALSE]}
    if(kind=='reverse')X<-X[,1,drop=FALSE]
    if(scaled){Y<-z(Y);X<-z(X);Z<-cbind(1,z(Z))}
    P<-diag(ncol(Z));if(scaled)P[1,1]<-0
    fit<-function(A)A-Z%*%solve(crossprod(Z)+ridge*P,crossprod(Z,A))
    U<-fit(Y);V<-fit(X);N<-nrow(U);d<-sum((crossprod(U,V)/N)^2)
    diagterm<-mean(rowSums(U^2)*rowSums(V^2))
    c(ordinary=d,diagonal_corrected=(N*d-diagterm)/(N-1))*(N/n)^2
  })
  list(values=vals,lambda=lam[-20])
}
cfg<-expand.grid(n=c(100,400,1600),D=c(7,21),scaled=c(FALSE,TRUE),kind=c('scalar','functional','reverse'),stringsAsFactors=FALSE)
results<-list()
if(Sys.getenv('SN_DIAG_REVERSE_ONLY')=='1')cfg<-cfg[cfg$kind=='reverse',] else cfg<-cfg[cfg$kind!='reverse',]
for(g in seq_len(nrow(cfg))) {
  a<-cfg[g,];delta<-pop(a$D,a$scaled,a$kind)
  sim<-parallel::mclapply(1:R,function(r) {
    dat<-sim_FARX1_diag(a$n,D=a$D,theta=.3,seed=20260918+g*100000+r)
    p<-path(dat,a$scaled,a$kind);v<-p$values;G<-v[,1:19,drop=FALSE]-v[,20]*rep(p$lambda^2,each=2)
    q<-sqrt(rowMeans(G^2));h<-apply(G,1,function(x)diff(range(x)))
    c(d_hat=v[1,20],d_corrected=v[2,20],quad=(v[1,20]-delta)/q[1]>cvtab[trim==.25 & method=='quad' & alpha==.05]$quantile,
      range=(v[1,20]-delta)/h[1]>cvtab[trim==.25 & method=='range' & alpha==.05]$quantile,
      quad_corrected=(v[2,20]-delta)/q[2]>cvtab[trim==.25 & method=='quad' & alpha==.05]$quantile,
      range_corrected=(v[2,20]-delta)/h[2]>cvtab[trim==.25 & method=='range' & alpha==.05]$quantile)
  },mc.cores=8)
  mat<-do.call(rbind,sim)
  if(any(!is.finite(mat)))stop('Failed diagnostic')
  row<-cbind(as.data.table(a),data.table(Delta=delta,R=R),as.data.table(as.list(colMeans(mat))))
  setnames(row,gsub('\\.ordinary|\\.diagonal_corrected','',names(row)))
  results[[g]]<-row;fwrite(rbindlist(results),file.path(out,if(Sys.getenv('SN_DIAG_REVERSE_ONLY')=='1')'reverse_diagnostics.csv' else 'finite_sample_diagnostics.csv'));print(row)
}
