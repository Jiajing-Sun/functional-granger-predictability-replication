#!/usr/bin/env Rscript
# Boundary calibration with standardized calendar HAR / reverse designs.
# Population thresholds use exact stationary lag covariances, not simulation.
source('FGC_simulation_functions_lambda20_plugin.R')
library(data.table)
cv<-fread('../processed/simulation/critical_values.csv');R<-1000L
pop<-function(error,kind) {
  psi<-if(error=='fMA1').7 else 0
  A<-rbind(c(.5,0,psi,0),c(.3,.4,0,psi),c(0,0,0,0),c(0,0,0,0))
  B<-rbind(c(1,0),c(0,1),c(1,0),c(0,1))
  P<-solve_discrete_lyapunov(A,B%*%t(B)/(1+psi^2))
  lag<-list(P);for(h in 1:22)lag[[h+1]]<-A%*%lag[[h]]
  cc<-function(i,li,j,lj)if(lj>=li)lag[[lj-li+1]][i,j] else lag[[li-lj+1]][j,i]
  if(kind=='HAR') {
    groups<-list(1,1:5,1:22)
    S<-outer(1:3,1:3,Vectorize(function(i,j)mean(outer(groups[[i]],groups[[j]],Vectorize(function(a,b)cc(2,a,2,b))))))
    cy<-vapply(groups,function(g)mean(vapply(g,function(a)cc(2,0,2,a),0.)),0.)
    cx<-vapply(groups,function(g)mean(vapply(g,function(a)cc(1,1,2,a),0.)),0.)
    C<-cc(2,0,1,1)-sum(cy*solve(S,cx))
  } else C<-cc(2,0,1,1)-cc(2,0,2,1)*P[2,1]/P[2,2]
  C^2/(P[1,1]*P[2,2])
}
z<-function(M)sweep(sweep(M,2,colMeans(M),'-'),2,apply(M,2,sd),'/')
configs<-list();j<-1L
add<-function(error='iid',kind='HAR',D=21,n=400,trim=.25,ridge=1,availability=1) {
  configs[[length(configs)+1L]]<<-data.table(error,kind,D,n,trim,ridge,availability)
}
for(e in c('iid','fMA1','t5'))for(k in c('HAR','reverse'))add(error=e,kind=k)
for(k in c('HAR','reverse')) {
  add(kind=k,D=7);add(kind=k,ridge=.1);add(kind=k,ridge=10)
  add(kind=k,trim=.5);add(kind=k,availability=.85)
  add(kind=k,n=1600,availability=.85)
}
add(kind='HAR',trim=0)
cfg<-unique(rbindlist(configs));results<-list()
for(g in seq_len(nrow(cfg))) {
  a<-cfg[g];delta<-pop(a$error,a$kind);lam<-a$trim+(1-a$trim)*(1:19)/20
  qv<-cv[trim==a$trim & method=='quad' & alpha==.05]$quantile
  hv<-cv[trim==a$trim & method=='range' & alpha==.05]$quantile
  sim<-parallel::mclapply(1:R,function(rr) {
    d<-sim_FARX1_diag(a$n+22,D=a$D,theta=.3,error_type=a$error,seed=20260920+g*100000+rr)
    ti<-23:(a$n+22);M<-rbinom(a$n+22,1,a$availability)==1
    if(a$kind=='HAR') {
      Y<-d$Y[ti,1,drop=FALSE];X<-d$X[ti-1,,drop=FALSE]
      Z<-cbind(d$Y[ti-1,1],vapply(ti,function(t)mean(d$Y[(t-5):(t-1),1]),0.),vapply(ti,function(t)mean(d$Y[(t-22):(t-1),1]),0.))
      valid<-M[ti-1]
    } else {
      Y<-d$Y[ti,,drop=FALSE];X<-d$X[ti-1,1,drop=FALSE];Z<-d$Y[ti-1,,drop=FALSE]
      valid<-M[ti]&M[ti-1]
    }
    vals<-vapply(c(lam,1),function(l) {
      k<-floor(a$n*l);ix<-which(valid & seq_len(a$n)<=k);N<-length(ix)
      if(N<ncol(Z)+4)return(NA_real_)
      ZZ<-cbind(1,z(Z[ix,,drop=FALSE]));YY<-z(Y[ix,,drop=FALSE]);XX<-z(X[ix,,drop=FALSE]);pen<-diag(ncol(ZZ));pen[1,1]<-0
      fit<-function(A)A-ZZ%*%solve(crossprod(ZZ)+a$ridge*pen,crossprod(ZZ,A))
      U<-fit(YY);V<-fit(XX);dd<-sum((crossprod(U,V)/N)^2)
      ((N*dd-mean(rowSums(U^2)*rowSums(V^2)))/(N-1))*(k/a$n)^2
    },0.)
    if(anyNA(vals))return(c(d=NA,quad=NA,range=NA))
    G<-vals[1:19]-lam^2*vals[20]
    c(d=vals[20],quad=(vals[20]-delta)/sqrt(mean(G^2))>qv,range=(vals[20]-delta)/diff(range(G))>hv)
  },mc.cores=8)
  mat<-do.call(rbind,sim)
  row<-cbind(a,data.table(Delta=delta,R=R,valid_replications=sum(complete.cases(mat))),as.data.table(as.list(colMeans(mat,na.rm=TRUE))))
  results[[g]]<-row;fwrite(rbindlist(results),'../processed/diagnostics/matched_diagnostics.csv');print(row)
}
