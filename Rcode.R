# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

r.tvec=function(data, otherX=NULL, lag=1, nbeta=50, ngamma=50, trim=0.05, print.out=FALSE){
 n=nrow(data)
 dY=data[(lag+2):n,]-data[(lag+1):(n-1),]
 dX=matrix(0,n-lag-1,2*lag)
 for (i in 1:lag){
  dX[,(2*i-1):(2*i)]=data[(lag+2-i):(n-i),]-data[(lag+1-i):(n-1-i),]
 }
 lm.coint=lm(data[,1]~data[,2]-1)
 beta.mean=lm.coint$coef[1]
 beta.sd=sqrt(diag(summary(lm.coint)$sigma*summary(lm.coint)$cov))[1]
 betagrid=matrix(seq(beta.mean-2*beta.sd,beta.mean+2*beta.sd,len=nbeta),ncol=1)
 X.1=embed(data,lag+2)[,3:4]
 ECT.1=X.1%*%c(1,-beta.mean)
 gamma.list=sort(unique(ECT.1))
 gammagrid=unique(gamma.list[round(seq(trim,1-trim,len=ngamma)*length(gamma.list))])
 beta.hat=NULL
 gamma.hat=NULL
 min.RSS=Inf
 for (i in 1:length(betagrid))
  for (j in 1:length(gammagrid)){
   ECTi=data[(lag+1):(n-1),1]-betagrid[i]*data[(lag+1):(n-1),2]
   d=ECTi<=gammagrid[j] 
   Z=cbind(ECTi*d,d,dX*d,otherX[(lag+1):(n-1),]*d,
ECTi*(1-d),1-d,dX*(1-d),otherX[(lag+1):(n-1),]*(1-d))
   A.hat=NULL
   A.hat=try(solve(t(Z)%*%Z)%*%t(Z)%*%dY,silent=TRUE)
   if ((is.numeric(A.hat)==TRUE)&&(sum(is.na(A.hat))==0)){
    RSS=sum((dY-Z%*%A.hat)^2)
    if (RSS<min.RSS){
     min.RSS=RSS
     beta.hat=betagrid[i]
     gamma.hat=gammagrid[j]
    }
   }
   else{
    RSS=Inf
   }   
   if (print.out==TRUE){
    cat("beta=",betagrid[i]," gamma=",gammagrid[j]," RSS=",RSS,"\n")
   }
  }
 ECTi=data[(lag+1):(n-1),1]-beta.hat*data[(lag+1):(n-1),2]
 d=ECTi<=gamma.hat 
 Z=cbind(ECTi*d,d,dX*d,otherX[(lag+1):(n-1),]*d,
ECTi*(1-d),1-d,dX*(1-d),otherX[(lag+1):(n-1),]*(1-d))
 A.hat=try(solve(t(Z)%*%Z)%*%t(Z)%*%dY,silent=TRUE)
 if (print.out==TRUE){
  cat("beta=",beta.hat,"\n")
  cat("gamma=",gamma.hat,"\n")
  cat("A=",A.hat,"\n")
 }
 list(beta=beta.hat,gamma=gamma.hat,A=A.hat,lag=lag,res=dY-Z%*%A.hat)
}

c.tvec=function(data, otherX=NULL, lag=1, nbeta=50, ngamma=50, trim=0.05){
 n=nrow(data)
 l=as.integer(lag)
 nb=as.integer(nbeta)
 ng=as.integer(ngamma)
 trm=as.double(trim)
 if (is.null(otherX))
  cX=as.integer(0)
 else
  cX=as.integer(ncol(otherX))
 d=data[,1:2]
 out=.C("tvec", as.double(d), as.double(otherX), n, l, nb, ng, trm, cX, beta=double(1), gamma=double(1), A=double(2*(4+4*l+2*cX)), res=double(2*(n-l-1)))
 list(beta=out$beta, gamma=out$gamma, A=matrix(out$A,ncol=2), lag=lag, res=matrix(out$res,ncol=2))
}

predict.tvec=function(tvec.model, last, dlast, other.last.X=NULL){
 lag=tvec.model$lag
 nox=length(other.last.X)
 w=last[1]-tvec.model$beta*last[2]
 gamma=tvec.model$gamma
 X=c(w,1)
 for (i in 1:lag)
  X=c(X,dlast[i,1],dlast[i,2])
 X=c(X,other.last.X)
 predicted=matrix(0,2,1)
 if (w<gamma)  
  predicted=t(tvec.model$A[1:(2*lag+2+nox),])%*%X
 else
  predicted=t(tvec.model$A[(2*lag+3+nox):(4*lag+4+2*nox),])%*%X
 predicted
}

r.tveci=function(data, lag=1, print.out=FALSE, icount.lim=5, count.lim=100, ep=1e-4, otherX=NULL,...){
 z=which(data[,3]==0)
 zz=which(data[,3]==1)
 p=data[,1:2]
 n=nrow(p)
  if (z[1]==1)
   p[1,2]=p[which(data[,3]!=0)[1],2]
  else{
   p[z[1],2]=p[z[1]-1,2]
  }
  if (length(z)>1){
   for (i in 2:length(z)){
    p[z[i],2]=p[z[i]-1,2]
   }
  }
 z=z[z>lag+1]
 nz=length(z)
 otv=r.tvec(p,lag=lag,...)
 idiff=Inf 
 idiffbest=Inf
 count=0
 icount=0
 pbest=p
 otvbest=otv
 while ((icount<=icount.lim)&(count<=count.lim)){
  old.p=p
  if (nz>0){
   for (i in 1:nz){
    zi=z[i]
    p[zi,2]=predict.tvec(otv,old.p[zi-1,],
matrix(old.p[(zi-1):(zi-lag),]-old.p[(zi-2):(zi-lag-1),],ncol=2))[2,1]+old.p[zi-1,2]
   }
  }
  idiff=sum(otv$res[,1]^2)
  otv=r.tvec(p,lag=lag,otherX=otherX,...)
  if (idiff < idiffbest-ep){
   icount=0
   idiffbest=idiff 
   pbest=p
   otvbest=otv  
  }
  else{
   icount=icount+1
  }
  count=count+1
 }
 if (count.lim<0)
  otvbest=r.tvec(p,lag=lag,otherX=otherX,...)
 list(beta=otvbest$beta, gamma=otvbest$gamma, A=otvbest$A, lag=lag, res=otvbest$res, imputed=pbest, otherX=otherX, steps=count)
}

c.tveci=function(data, lag=1, print.out=FALSE, icount.lim=5, count.lim=100, ep=1e-4, otherX=NULL, ...){
 z=which(data[,3]==0)
 zz=which(data[,3]==1)
 p=data[,1:2]
 n=nrow(p)
  if (z[1]==1)
   p[1,2]=p[which(data[,3]!=0)[1],2]
  else{
   p[z[1],2]=p[z[1]-1,2]
  }
  if (length(z)>1){
   for (i in 2:length(z)){
    p[z[i],2]=p[z[i]-1,2]
   }
  }
 z=z[z>lag+1]
 nz=length(z)
 otv=c.tvec(p,lag=lag,otherX=otherX,...)
 idiff=Inf 
 idiffbest=Inf
 count=0
 icount=0
 pbest=p
 otvbest=otv
 while ((icount<=icount.lim)&(count<=count.lim)){
  old.p=p
  if (nz>0){
    for (i in 1:nz){
     zi=z[i]
     p[zi,2]=predict.tvec(otv,old.p[zi-1,],
matrix(old.p[(zi-1):(zi-lag),]-old.p[(zi-2):(zi-lag-1),],ncol=2))[2,1]+old.p[zi-1,2]
    }
  }
  idiff=sum(otv$res[,1]^2)
  otv=c.tvec(p,lag=lag,otherX=otherX,...)
  if (idiff < idiffbest-ep){
   icount=0
   idiffbest=idiff 
   pbest=p
   otvbest=otv  
  }
  else{
   icount=icount+1
  }
  count=count+1
 }
 if (count.lim<0)
  otvbest=c.tvec(p,lag=lag,...)
 list(beta=otvbest$beta, gamma=otvbest$gamma, A=otvbest$A, lag=lag, res=otvbest$res, imputed=pbest, otherX=otherX, steps=count)
}

predict.tveci=function(tveci.model){
 nox=length(tveci.model$otherX)
 n=nrow(tveci.model$imputed)
 w=tveci.model$imputed[n,1]-tveci.model$beta*tveci.model$imputed[n,2]
 gamma=tveci.model$gamma
 X=c(w,1)
 dlast=matrix(tveci.model$imputed[n:(n-tveci.model$lag+1),],tveci.model$lag,2)-matrix(tveci.model$imputed[(n-1):(n-tveci.model$lag),],tveci.model$lag,2)
 for (i in 1:tveci.model$lag)
  X=c(X,dlast[i,1],dlast[i,2])
 X=c(X,tveci.model$otherX[n,])
 predicted=matrix(0,2,1)
 if (w<tveci.model$gamma)  
  predicted=t(tveci.model$A[1:(2*tveci.model$lag+2+nox),])%*%X
 else
  predicted=t(tveci.model$A[(2*tveci.model$lag+3+nox):(4*tveci.model$lag+4+2*nox),])%*%X
 predicted
}

