generate_data <- function(l, path_name) {
  set.seed(100+l)
  Z<-matrix(rbinom(n,size=1,prob=0.5),ncol=1)
  dat<-sim.data_dep(n,Z,a,b,c,Delta,Xi1,Sigma2)
  M<-dat$M
  R<-dat$R
  
  colnames(M)<-paste0("M",1:k)
  dd0<-data.frame(Z=Z,M,R=R)
  
  # standardize data
  m.Z<-mean(Z)
  m.M<-apply(M,2,mean)
  m.R<-mean(R)
  sd.Z<-sd(Z)
  sd.M<-apply(M,2,sd)
  sd.R<-sd(R)
  
  Z<-scale(Z)
  M<-scale(M)
  R<-scale(R)
  dd<-data.frame(Z=Z,M,R=R)
  
  save(list=c("dd0","dd","dat","Z","M","R","m.Z","m.M","m.R","sd.Z","sd.M","sd.R"),file = path_name)
}
runonce<-function(l) { # l could be the number of replicates and is used to decide the seed
  # Generate data
  if(file.exists(paste0(dir.data,"/RUN_",l,"_Data.RData"))==FALSE) { generate_data(l, paste0(dir.data,"/RUN_",l,"_Data.RData")) }
  load(paste0(dir.data,"/RUN_",l,"_Data.RData"))
  
  run.file<-paste0(dir.data,"/RUN_",l,"_",method,".RData")
  
  if(file.exists(run.file)==FALSE) {
    # calculate scale back multiplier
    A.pf<-sd.Z/sd.M
    B.pf<-sd.M/sd.R
    C.pf<-sd.Z/sd.R
    AB.pf<-A.pf*B.pf
    
    re<-vector("list",length=length(lambda))
    AB.es=A.estt=B.est<-matrix(NA,k,length(lambda))
    C.est<-rep(NA,length(lambda))
    
    if(method=="Lasso") {
      for(i in length(lambda):1) {
        out<-NULL
        if(i==length(lambda)) {
          # starting from the largest lambda value
          try(out<-mediation_net_ADMM_NC(Z,M,R,
                                         lambda=0,
                                         omega=lambda[i],
                                         phi=phi,
                                         Phi1=Phi1,Phi2=Phi2,
                                         rho=rho,rho.increase=FALSE,
                                         tol=tol,max.itr=max.itr,
                                         thred=thred,
                                         Sigma1=Sigma10,Sigma2=Sigma20,trace=FALSE))
        } else {
          # for smaller lambda (ith lambda), use the (i+1)th lambda results as burn-in 
          try(out<-mediation_net_ADMM_NC(Z,M,R,lambda=0,omega=lambda[i],phi=phi,Phi1=Phi1,Phi2=Phi2,
                                         rho=rho,rho.increase=FALSE,tol=tol,max.itr=max.itr,thred=thred,Sigma1=Sigma10,Sigma2=Sigma20,trace=FALSE,
                                         Theta0=matrix(c(1,A.est[,i+1]*(sd.Z/sd.M)),nrow=1),D0=matrix(c(C.est[i+1]*(sd.Z/sd.R),B.est[,i+1]*(sd.M/sd.R)),ncol=1),
                                         alpha0=matrix(c(1,A.est[,i+1]*(sd.Z/sd.M)),nrow=1),beta0=matrix(c(C.est[i+1]*(sd.Z/sd.R),B.est[,i+1]*(sd.M/sd.R)),ncol=1)))
        }
        
        if(is.null(out)==FALSE) {
          re[[i]]<-out
          # scale the estimate back to the original scale
          B.est[,i]<-out$B*(sd.R/sd.M)
          C.est[i]<-out$C*(sd.R/sd.Z)
          A.est[,i]<-out$A*(sd.M/sd.Z)
          AB.est[,i]<-A.est[,i]*B.est[,i]
        }
        
        print(paste0("lambda index ",i))
      }
    }
    if(length(grep("PathLasso-",method))>0) {
      omega.p.idx<-as.numeric(sub("PathLasso-","",method))
      
      for(i in length(lambda):1) {
        out<-NULL
        if(i==length(lambda)) {
          # starting from the largest lambda value
          try(out<-mediation_net_ADMM_NC(Z,M,R,lambda=lambda[i],omega=omega.p[omega.p.idx]*lambda[i],phi=phi,Phi1=Phi1,Phi2=Phi2,
                                         rho=rho,rho.increase=FALSE,tol=tol,max.itr=max.itr,thred=thred,Sigma1=Sigma10,Sigma2=Sigma20,trace=FALSE))
        } else {
          # for smaller lambda (ith lambda), use the (i+1)th lambda results as burn-in 
          try(out<-mediation_net_ADMM_NC(Z,M,R,lambda=lambda[i],omega=omega.p[omega.p.idx]*lambda[i],phi=phi,Phi1=Phi1,Phi2=Phi2,
                                         rho=rho,rho.increase=FALSE,tol=tol,max.itr=max.itr,thred=thred,Sigma1=Sigma10,Sigma2=Sigma20,trace=FALSE,
                                         Theta0=matrix(c(1,A.est[,i+1]*(sd.Z/sd.M)),nrow=1),D0=matrix(c(C.est[i+1]*(sd.Z/sd.R),B.est[,i+1]*(sd.M/sd.R)),ncol=1),
                                         alpha0=matrix(c(1,A.est[,i+1]*(sd.Z/sd.M)),nrow=1),beta0=matrix(c(C.est[i+1]*(sd.Z/sd.R),B.est[,i+1]*(sd.M/sd.R)),ncol=1)))
        }
        
        if(is.null(out)==FALSE)
        {
          re[[i]]<-out
          
          # scale the estimate back to the original scale
          B.est[,i]<-out$B*(sd.R/sd.M)
          C.est[i]<-out$C*(sd.R/sd.Z)
          A.est[,i]<-out$A*(sd.M/sd.Z)
          AB.est[,i]<-A.est[,i]*B.est[,i]
        }
        
        print(paste0("lambda index ",i))
      }
    }
    
    save(list=c("re","AB.est","A.est","B.est","C.est","A.pf","B.pf","C.pf","AB.pf"),file=run.file)
  }
}

runonce.KSC<-function(l) {
  if(!file.exists(paste0(dir.data,"/RUN_",l,"_Data.RData"))) { generate_data(l, paste0(dir.data,"/RUN_",l,"_Data.RData")) }
  load(paste0(dir.data,"/RUN_",l,"_Data.RData"))
 
  run.file<-paste0(dir.data,"/RUN_",l,"_KSC_",method,".RData")
  
  if(file.exists(run.file)==FALSE) {
    if(method=="Lasso") {
      out<-NULL
      try(out<-mediation_net_ADMM_NC_KSC(Z,M,R,zero.cutoff=zero.cutoff,n.rep=n.rep,vss.cut=vss.cut,lambda=0,omega=lambda,
                                         phi=phi,Phi1=Phi1,Phi2=Phi2,rho=rho,rho.increase=FALSE,tol=tol,max.itr=max.itr,thred=thred,
                                         Sigma1=Sigma10,Sigma2=Sigma20,trace=FALSE,Theta0=NULL,D0=NULL,alpha0=NULL,beta0=NULL))
    }
    if(length(grep("PathLasso-",method))>0) {
      omega.p.idx<-as.numeric(sub("PathLasso-","",method))
      out<-NULL
      try(out<-mediation_net_ADMM_NC_KSC(Z,M,R,zero.cutoff=zero.cutoff,n.rep=n.rep,vss.cut=vss.cut,lambda=lambda,omega=omega.p[omega.p.idx]*lambda,
                                         phi=phi,Phi1=Phi1,Phi2=Phi2,rho=rho,rho.increase=FALSE,tol=tol,max.itr=max.itr,thred=thred,
                                         Sigma1=Sigma10,Sigma2=Sigma20,trace=FALSE,Theta0=NULL,D0=NULL,alpha0=NULL,beta0=NULL))
    } 
    save(list=c("out"),file=run.file)
  }
}