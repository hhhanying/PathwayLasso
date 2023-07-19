##########################################################
# Pathway lasso
##########################################################

rm(list=ls())

source("functions.R")
source("ADMM_adp_functions.R")
source("VSS_adp_functions.R")

##########################################################
# parameter setting
k<-50

n<-50

p<-1

sigma20_1<-200^2
sigma20_2<-200^2

A.nz<-c(3,9,18,36,0,0,12,4)*c(c(1,-1,1,-1),c(1,1,-1,1))
B.nz<-c(12,4,2,1,4,9,0,0)*c(c(1,1,-1,1),c(1,-1,-1,1))
A<-matrix(c(A.nz,rep(0,k-length(A.nz))),nrow=1)
B<-matrix(c(B.nz,rep(0,k-length(B.nz))),ncol=1)
AB<-t(A)*B
C<-max(abs(AB))
C2<-C+A%*%B
D<-rbind(C,B)
##########################################################

##########################################################
# Pathway Lasso method parameters
phi<-2

rho<-1             # ADMM parameter
max.itr<-5000
tol<-1e-6
thred<-1e-6

thred2<-1e-3

omega.p<-c(0,0.1,1)      # omega = omega.p*lambda

# tuning parameter lambda
lambda<-c(10^c(seq(-5,-3,length.out=5),seq(-3,0,length.out=21)[-1],seq(0,2,length.out=11)[-1],seq(2,4,length.out=6)[-1])) # lambda values

# tuning parameter selection by variable selection stability
# kappa selection criterion parameter
n.rep<-5
vss.cut<-0.1
##########################################################

L<-200

dir.data0<-paste0(getwd(),"/phi_",phi)
if(file.exists(dir.data0)==FALSE) { dir.create(dir.data0) }

method0<-c("Lasso",paste0("PathLasso-",1:length(omega.p)))

# add dependence between mediators
rho.M<-c(-0.4,0,0.4)

for(mm in 1:length(rho.M)) {
  # Covariance matrix
  Sigma1<-matrix(0,k,k)
  set.seed(10000)
  Sigma1[upper.tri(Sigma1,diag=FALSE)]<-rho.M[mm]*sigma20_1*rbinom(length(Sigma1[upper.tri(Sigma1,diag=FALSE)]),size=1,prob=1/k)
  Sigma1<-Sigma1+t(Sigma1)
  diag(Sigma1)<-rep(sigma20_1,k)
  Sigma2<-sigma20_2
  
  Sigma<-cbind(rbind(Sigma1,rep(0,k)),c(rep(0,k),Sigma2))
  
  Sigma1.chol<-chol(Sigma1)
  Xi1<-diag(diag(Sigma1.chol)^2)
  Delta<-diag(rep(1,k))-solve(Sigma1.chol/diag(Sigma1.chol))
  
  # 
  Sigma10<-diag(rep(1,k))
  Sigma20<-matrix(1,1,1)
  
  a<-A%*%(diag(rep(1,k))-Delta)
  b<-B
  c<-C
  
  tau2<-(t(B)%*%Sigma1%*%B+Sigma2)[1,1]
  Rho<-Sigma1%*%B
  
  fname<-paste0("Sim_k",k,"_corM",rho.M[mm])
  dir.data<-paste0(dir.data0,"/",fname)
  if(file.exists(dir.data)==FALSE)
  {
    dir.create(dir.data)
  }
  
  require(parallel)
  
  for(ii in 1:length(method0)) {
    method<-method0[ii]    
    mclapply(1:L,runonce,mc.cores=20)
  }
  for(ii in 1:length(method0)) {
    method<-method0[ii]
    mclapply(1:L,runonce.KSC,mc.cores=20)
  }
}


