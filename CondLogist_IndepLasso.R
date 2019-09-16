

CondLogist_IndepLasso <- function(X,y,pairs,strata,method=c("BIC","BIC-R","CV","CV-OSL")){
  NumFolds=10
  W=NULL
  Criteria=method
  Test=length(which(method %in% c("BIC","BIC-R","CV","CV-OSL")))
  if (Test==0) {
    cat("Please choose at least one of the following criteria: BIC, BIC-R, CV and CV-OSL ")  
  }else{
    Y=y
    X=X
    p=ncol(X)
    N=length(Y)
    K=max(strata)
    ncs=table(strata)
    ORD=order(strata,pairs,-y)
    
  
  Y1=Y[ORD]
  Y2=Y1
  Y2[which(Y2!=0)]=1
  X1=X[ORD,]
  Z1=NULL
  for (i in 1:K) {
    Z1=c(Z1,rep(i,ncs[i]))
  }


  if(is.null(W)){W=matrix(rep(1,p*K),ncol=K)}


  ALL_COEFS=ALL_COEFSIter=array(NA,c(p,50,K))
  ALL_COEFSR=list(NULL)
  ALL_CVMEAN=ALL_CVSE=ALL_CVMEANIter=ALL_CVSEIter=ALL_BIC=ALL_BICR=matrix(NA,ncol = 50,nrow = K)
  bic=cvmean=cv_1se=cv_1.39se=cv_1.96se=bicR=cvmeaniter=cv_1se_iter=cv_1.39se_iter=cv_1.96se_iter=NULL
  verif=verifI=list()
  for (k in 1:K) {
    cat(k)
    y_k=Y2[which(Z1==k)]
    x_k=X1[which(Z1==k),]
    BIC=NULL
    paire=seq(2,length(y_k),by = 2)
    impaire=seq(1,length(y_k),by = 2)
    
    MatDiagWeightsInv = sparseMatrix(i = 1:(p), j = 1:(p) ,x= 1/c(W[,k]))
    
    x_k=x_k%*%MatDiagWeightsInv  
    
    
    mod0         =  clogitLasso(x_k,y_k,strata = rep(1:(length(y_k)/2),each=2),nbfraction = 50,epsilon = 0.001)
    ALL_COEFS[,,k]=t(mod0$beta)
    if ("CV" %in% Criteria || "CV-OSL" %in% Criteria) {
      
    mod1         =  CV.clogit(cbind(y_k,x_k,y_k),NumFolds,mod0$fraction)
    
    eval(parse(text = paste0("verif$k",k,"=mod1$verif")))
    
    ALL_CVMEAN[k,]=mod1$meancv
    ALL_CVSE[k,]=mod1$secv
    }
    
    if ("BIC" %in% Criteria) {
    for (l in 1:dim(mod0$beta)[1]) {
      XB1=exp(as.vector(x_k[impaire,]%*%mod0$beta[l,]))
      XB2=exp(as.vector(x_k[paire,]%*%mod0$beta[l,]))
      LX=XB1/(XB1+XB2)
      LL = sum(log(LX))
      LL=as.numeric(LL)
      
      BICC= -2*LL + length(which(mod0$beta[l,]!=0))*log(length(y_k))
      
      BIC=c(BIC,BICC)
    }
      ALL_BIC[k,]=BIC
    
    bic=cbind(bic,as.vector(MatDiagWeightsInv%*%as.vector(mod0$beta[which.min(BIC),])))
    }
    if ("CV" %in% Criteria || "CV-OSL" %in% Criteria) {
    cvmean=cbind(cvmean,as.vector(MatDiagWeightsInv%*%as.vector(mod0$beta[which.min(mod1$meancv),])))
    
    for (cc in c(1,1.39,1.96)) {
   
    Cond= mod1$meancv[which.min(mod1$meancv)] + mod1$secv[which.min(mod1$meancv)]*cc
    POsition=which(mod1$meancv<=Cond)
    if (length(POsition) >1) {
    Comp=apply(mod0$beta[POsition,], 1, function(b){
      length(which(b!=0))
    })
    
    OneSe=POsition[which.min(Comp)]
    }else{OneSe=POsition}
    
    

    eval(parse(text = paste0("cv_",cc,"se  =cbind(cv_",cc,"se,as.vector(MatDiagWeightsInv%*%as.vector(mod0$beta[OneSe,])))")))
    }
    }
    
    if ("BIC-R" %in% Criteria) {
    relax = sel.BIC.relaxed_IND(t(mod0$beta), x_k, y_k ,rep(1:(length(y_k)/2),each=2), p)
    bicR=cbind(bicR,as.vector(MatDiagWeightsInv%*%as.vector(relax[[1]][[which.min(relax[[2]])]]$coef.relax)))
    eval(parse(text = paste0("ALL_COEFSR$Stratum_",k,"=relax[[1]]")))
    ALL_BICR[k,]=relax[[2]]
    }
    
    
    
    
  }
  
  if ("CV-OSL" %in% Criteria) {
  
  for (k in 1:K) {
    cat(k)
    y_k=Y2[which(Z1==k)]
    x_k=X1[which(Z1==k),]


    paire=seq(2,length(y_k),by = 2)
    impaire=seq(1,length(y_k),by = 2)
    
    MatDiagWeightsInv = sparseMatrix(i = 1:(p), j = 1:(p) ,x= 1/c(W[,k]))
    if (K==1) {
      WIter= 1/(abs(c(cvmean))+1e-4) 
    }else{
      WIter= 1/(abs(c(cvmean[,k]))+1e-4)}
    MatDiagWeightsIterInv = sparseMatrix(i = 1:(p), j = 1:(p) ,x= 1/WIter)
    
    x_kk=x_k%*%MatDiagWeightsInv  
    x_kk=x_kk%*%MatDiagWeightsIterInv
    
   
    mod0         =  clogitLasso(x_kk,y_k,strata = rep(1:(length(y_k)/2),each=2),nbfraction = 50,epsilon = 0.001)
    ALL_COEFSIter[,,k]=t(mod0$beta)

    
    mod1         =  CV.clogit1(cbind(y_k,x_k%*%MatDiagWeightsInv,y_k),NumFolds,mod0$fraction)
    eval(parse(text = paste0("verifI$k",k,"=mod1$verif")))
    CVMEANIter= mod1$meancv
    CVSEIter=mod1$secv
      
    ALL_CVMEANIter[k,]=CVMEANIter
    ALL_CVSEIter[k,]=CVSEIter
    
    
    cvmeaniter=cbind(cvmeaniter,as.vector(MatDiagWeightsIterInv%*%MatDiagWeightsInv%*%as.vector(mod0$beta[which.min(CVMEANIter),])))
    
    for (cc in c(1,1.39,1.96)) {
      
    Cond= CVMEANIter[which.min(CVMEANIter)] + CVSEIter[which.min(CVMEANIter)]*cc
    POsition=which(CVMEANIter<=Cond)
    if (length(POsition)>1) {
        
    
    Comp=apply(mod0$beta[POsition,], 1, function(b){
      length(which(b!=0))
    })
    OneSeiter=POsition[which.min(Comp)]
    }else{OneSeiter=POsition}
    
    eval(parse(text = paste0("cv_",cc,"se_iter  =cbind(cv_",cc,"se_iter,as.vector(MatDiagWeightsIterInv%*%MatDiagWeightsInv%*%as.vector(mod0$beta[OneSeiter,])))")))
    }
    
    
    
  }
  }
  bic=bic;bicR=bicR;cvmean=cvmean;cv_1se=cv_1se;cv_1.39se=cv_1.39se;cv_1.96se=cv_1.96se;cvmeaniter=cvmeaniter;cv_1se_iter=cv_1se_iter;cv_1.39se_iter=cv_1.39se_iter;cv_1.96se_iter=cv_1.96se_iter;
       ALL_COEFS=ALL_COEFS;ALL_BIC=ALL_BIC;ALL_CVMEAN=ALL_CVMEAN;ALL_CVSE=ALL_CVSE;ALL_COEFSR=ALL_COEFSR;
       ALL_BICR=ALL_BICR;ALL_COEFSIter=ALL_COEFSIter;ALL_CVMEANIter=ALL_CVMEANIter;ALL_CVSEIter=ALL_CVSEIter;
  
  return(list(BIC=bic,BIC_R=bicR,CV=cvmean,CV.1se=cv_1se,CV_OSL=cvmeaniter,CV.1se_OSL=cv_1se_iter))
  }
}
