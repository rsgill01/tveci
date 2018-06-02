// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.

#include <R.h>

extern int dgels_(char *trans, int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *work, int *lwork, int *info);

extern int dcopy_(int *n, double *dx, int *incx, double *dy, int *incy);

extern int dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);

void merge(double *A,double *L,double *R,int leftCount,int rightCount){
 int i=0;
 int j=0;
 int k=0;
 while (i<leftCount && j<rightCount){
  if (L[i] < R[j]) 
   A[k++]=L[i++];
  else A[k++]=R[j++];
 }
 while (i<leftCount) 
  A[k++]=L[i++];
 while (j<rightCount) 
  A[k++]=R[j++];
}

void sort(double *A,int n){
 int mid,i;
 double *L;
 double *R;
 if (n < 2) 
  return;
 mid=n/2;
 L=(double*)malloc(mid*sizeof(double)); 
 R=(double*)malloc((n-mid)*sizeof(double)); 
 for (i=0;i<mid;i++) 
  L[i]=A[i];
 for (i=mid;i<n;i++) 
  R[i-mid]=A[i];
 sort(L,mid);
 sort(R,n-mid);
 merge(A,L,R,mid,n-mid); 
 free(L);
 free(R);
}

void unique(double *A,int *n){
 int i=1;
 int j=1;
 while (i<*n){
  if (A[j++]==A[i-1])
   (*n)--;
  else
   A[i++]=A[j-1];
 }
}

void tvec(double *d, double *otherX, int *n, int *l, int *nb, int *ng, double *trm, int *cX, double *beta, double *gamma, double *A, double *res){
 int i,j,k;
 int kni,kn,jni,uni,tempint;
 double *ECbeta;
 double *ECbetaSD;
 double *betagrid;
 double *gammagrid;
 double *dY;
 double *dYtemp;
 double *dX;
 double *dZ;
 double *dZtemp;
 double *Atemp;
 double *tempres;
 double *ECTi;
 double *ECTm1;
 double *temprY;
 double *temprX;
 double *wrk;
 double *wrkr;
 char trans='N';
 int intone=1;
 int inttwo=2;
 double done=1.0;
 double dmone=-1.0;
 int info;
 int ni=*n-*l-1;
 int twni=2*ni;
 int nj=4*(*l)+2**cX+4;
 int twnj=2*nj;
 int ninj=ni*nj;
 int lw=nj+160*ni;
 int lwr=1+16*(*n);
 int ib;
 int ig;
 double tempsum=0.0;
 double tempri;
 double minrss=-1.0;
 dY=malloc(2*ni*sizeof(double));
 dYtemp=malloc(2*ni*sizeof(double));
 dX=malloc(2*(*l)*ni*sizeof(double));
 dZ=malloc(ninj*sizeof(double));
 dZtemp=malloc(ninj*sizeof(double));
 Atemp=malloc(2*nj*sizeof(double));
 tempres=malloc(2*ni*sizeof(double));
 ECTi=malloc(ni*sizeof(double));
 ECTm1=malloc(2*(*l)*ni*sizeof(double));
 temprY=malloc((*n)*sizeof(double));
 temprX=malloc((*n)*sizeof(double));
 wrk=malloc(lw*sizeof(double));
 wrkr=malloc(lwr*sizeof(double));
 ECbeta=malloc(sizeof(double));
 ECbetaSD=malloc(sizeof(double));
 betagrid=malloc((*nb)*sizeof(double));
 gammagrid=malloc((*ng)*sizeof(double));

 for (k=0;k<2;k++){
  kni=k*ni;
  kn=k*(*n);
  for (i=0;i<ni;i++)
   dY[i+kni]=d[kn+i+*l+1]-d[kn+i+*l];
  for (j=0;j<*l;j++){
   jni=2*j*ni;
   for (i=0;i<ni;i++)
    dX[i+kni+jni]=d[kn+i+*l-j]-d[kn+i+*l-j-1];
  }
 }
 dcopy_(n,d,&intone,temprY,&intone);
 for (i=0;i<*n;i++)
  temprX[i]=d[*n+i];
 dgels_(&trans,n,&intone,&intone,temprX,n,temprY,n,wrkr,&lwr,&info);
 *ECbeta=temprY[0];
 for (i=0;i<*n;i++){
  tempri=d[i]-d[*n+i]*(*ECbeta);
  tempsum+=tempri*tempri;
 }
 *ECbetaSD=sqrt(tempsum/(*n-1));
 tempsum=0.0;
 for (i=0;i<*n;i++)
  tempsum+=d[*n+i]*d[*n+i];
 *ECbetaSD=sqrt(*ECbetaSD/tempsum);
 double betagridstart=*ECbeta-2*(*ECbetaSD);
 double betatick=4*(*ECbetaSD)/(*nb-1);
 for (i=0;i<*nb;i++)
  betagrid[i]=betagridstart+betatick*i;
 for (i=0;i<ni;i++)
  ECTm1[i]=d[i+*l]-(*ECbeta)*d[i+*l+*n];
 sort(ECTm1,ni);
 uni=ni;
 unique(ECTm1,&uni);
 double gammatick=(1-2*(*trm))/(*ng-1);
 for (i=0;i<*ng;i++){
  tempint=(int)(round(uni*(*trm+gammatick*i)));
  tempint--;
  if (tempint>=0)
   gammagrid[i]=ECTm1[tempint];
 }
 unique(gammagrid,ng);
 for (ib=0;ib<*nb;ib++){
  for (i=0;i<ni;i++)
   ECTi[i]=d[i+*l]-d[*n+i+*l]*betagrid[ib];	   
  for (ig=0;ig<*ng;ig++){
   for (i=0;i<ni;i++){
    if (ECTi[i]<=gammagrid[ig]){
     dZ[i]=ECTi[i];
     dZ[i+ni]=1;
     for (j=0;j<*l;j++){
      dZ[i+ni*(2+2*j)]=dX[i+ni*2*j];
      dZ[i+ni*(3+2*j)]=dX[i+ni*(2*j+1)];
     }      
     for (j=0;j<*cX;j++)
      dZ[i+ni*(2+2**l+j)]=otherX[(*l)+i+j*(*n-1)];            
     dZ[i+ni*(2**l+*cX+2)]=0;
     dZ[i+ni*(2**l+*cX+3)]=0;
     for (j=0;j<*l;j++){
      dZ[i+ni*(2**l+*cX+4+2*j)]=0;
      dZ[i+ni*(2**l+*cX+5+2*j)]=0;
     }      
     for (j=0;j<*cX;j++)
      dZ[i+ni*(4+4**l+*cX+j)]=0;            
    }
    else{
     dZ[i]=0;
     dZ[i+ni]=0;
     for (j=0;j<*l;j++){
      dZ[i+ni*(2+2*j)]=0;
      dZ[i+ni*(3+2*j)]=0;
     }      
     for (j=0;j<*cX;j++)
      dZ[i+ni*(2+2**l+j)]=0;  
     dZ[i+ni*(2**l+*cX+2)]=ECTi[i];
     dZ[i+ni*(2**l+*cX+3)]=1;
     for (j=0;j<*l;j++){
      dZ[i+ni*(2**l+*cX+4+2*j)]=dX[i+ni*2*j];
      dZ[i+ni*(2**l+*cX+5+2*j)]=dX[i+ni*(2*j+1)];
     }      
     for (j=0;j<*cX;j++)
      dZ[i+ni*(4+4**l+*cX+j)]=otherX[(*l)+i+j*(*n-1)];
    }
   }
   dcopy_(&twni,dY,&intone,dYtemp,&intone);
   dcopy_(&ninj,dZ,&intone,dZtemp,&intone);
   dgels_(&trans,&ni,&nj,&inttwo,dZtemp,&ni,dYtemp,&ni,wrk,&lw,&info);
   if (info==0){
    dcopy_(&nj,dYtemp,&intone,Atemp,&intone);
    dcopy_(&nj,dYtemp+ni,&intone,Atemp+nj,&intone);
    dcopy_(&twni,dY,&intone,dYtemp,&intone);
    dgemm_(&trans,&trans,&ni,&inttwo,&nj,&dmone,dZ,&ni,Atemp,&nj,&done,dYtemp,&ni);
    dcopy_(&twni,dYtemp,&intone,tempres,&intone);
    tempsum=0.0;
    for (i=0;i<ni;i++)
     tempsum+=tempres[i]*tempres[i]+tempres[i+ni]*tempres[i+ni];
    if ((minrss<0)|(tempsum<minrss)){
     minrss=tempsum;
     *beta=betagrid[ib];
     *gamma=gammagrid[ig];
     dcopy_(&twnj,Atemp,&intone,A,&intone);
     dcopy_(&twni,tempres,&intone,res,&intone);
    } 
   }
  }
 }
 free(gammagrid);
 free(betagrid);
 free(ECbetaSD);
 free(ECbeta);
 free(wrkr);
 free(wrk);
 free(tempres);
 free(Atemp);
 free(temprX);
 free(temprY);
 free(ECTm1);
 free(ECTi);
 free(dZtemp);
 free(dZ);
 free(dX);
 free(dYtemp);
 free(dY);
}
