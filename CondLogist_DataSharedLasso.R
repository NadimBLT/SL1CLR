Calise = function(Feat, ncs)
{
  # block diagonal bg(X1,X2,....,XK)
  nstr = length(ncs)
  n    = nrow(Feat)
  p    = ncol(Feat)
  Fin  = matrix(0, n, p*nstr)
  i0   = 1
  for (c in 1:nstr)
  {
    Fin[i0:(i0-1+ncs[c]), ((c-1)*p+1):(c*p)] = Feat[i0:(i0-1+ncs[c]),]
    i0 										 = i0 + ncs[c]
  }
  return(Fin)
}


Hierise = function(Feat, ncs)
{
  # block diagonal bg(X1,X2,....,XK) puis 1{C=1},...,1{C=K} puis X, puis 1
  Fin = cbind(Feat, Calise(Feat,ncs))
  return(Fin)
}

mise_forme_coef_hiers = function(coefs, nstr)
{
  p        = length(coefs)/(nstr+1) 
  Coefbar  = coefs[1:p]
  Coeforme = NULL
  coefsred = coefs[-(1:p)]
  for (c in 1:nstr)
  {
    Coeforme = c(Coeforme, Coefbar + coefsred[((c-1)*p+1):(c*p)])
  }
  return(Coeforme)	
}

refo <- function(vec,K){
  mat=matrix(vec, ncol = K+1,byrow = FALSE)
  mat1 = (mat[,-1] + matrix(rep(mat[,1],K),ncol = K,byrow = FALSE))
  return(mat1)
}


CondLogist_DataSharedLasso <- function(X,y,pairs,matchSet,method=c("BIC","BIC-R","CV","CV-OSL")){
  strata=matchSet
  WL1=NULL
  WL2=NULL
  SeqTau=1
  tauks=1
  balanced=TRUE
  NumFolds=10
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
  
  Xx=Hierise(X1,ncs)
  ratio=0.001
  Nlam=50
  
  if(tauks==1){tauks=rep(1,K)
  if(balanced== FALSE){tauks=sqrt(K)*sqrt(ncs/sum(ncs))}
  }
  
  if(is.null(WL1)){WL1=rep(1,p)}
  if(is.null(WL2)){WL2=rep(tauks,each=p)}


  MatDiagWeightsInv = sparseMatrix(i = 1:(p*(K+1)), j = 1:(p*(K+1)) ,x= 1/c(WL1,WL2))


  if(SeqTau!="SEQ50"){ratio_lam1_lam2=SeqTau}else{
  ratio.max       = 1.01*sum(tauks)
  ratio_lam1_lam2 = seq(ratio.max,ratio.max*ratio, length.out=Nlam)
  }
  
  ALL_COEFS=NULL
  BIC=CVMEAN=CVSE=NULL
  bic=cvmean=cv_1se=cv_1.39se=cv_1.96se=bicR=NULL

  paire=seq(2,N,by = 2)
  impaire=seq(1,N,by = 2)
  for (rk in ratio_lam1_lam2) {

    rk_mat = sparseMatrix(i = 1:(p*(K+1)), j = 1:(p*(K+1)) , x=c(rep(1/rk,p), rep(1,K*(p))))
    calX.Hier.rat = Xx%*%MatDiagWeightsInv%*%rk_mat


    mod0         =  clogitLasso(calX.Hier.rat,Y2,strata = pairs,nbfraction = 50,epsilon = ratio)
    
    if ("CV" %in% Criteria || "CV-OSL" %in% Criteria) {
      
    mod1        =  CV.clogit(cbind(Y2,calX.Hier.rat,Y1),NumFolds,mod0$fraction)
    verif=mod1$verif
}
    
 if ("BIC" %in% Criteria) {
   
    for (l in 1:dim(mod0$beta)[1]) {
    XB1=exp(as.vector(calX.Hier.rat[impaire,]%*%mod0$beta[l,]))
    XB2=exp(as.vector(calX.Hier.rat[paire,]%*%mod0$beta[l,]))
    LX=XB1/(XB1+XB2)
    LL = sum(log(LX))
    LL = as.numeric(LL)
    
    BICC= -2*LL + length(which(mod0$beta[l,]!=0))*log(length(Y2))
    
    BIC=c(BIC,BICC)
    }
   
   
}
    
    if ("CV" %in% Criteria || "CV-OSL" %in% Criteria) {
      
   
    CVMEAN= c(CVMEAN, mod1$meancv)
    
    CVSE=c(CVSE, mod1$secv)
    }
    
    bbeta		 = t(mod0$beta)
    
    
    ALL_COEFS = cbind(ALL_COEFS, apply(bbeta,2,function(v){as.vector(Matrix(diag(c(rep(1/rk,p), rep(1,K*(p)))), sparse=T)%*%MatDiagWeightsInv%*%v)}))
    
    

  }
  
  

  relaxed=list(NULL,
               NULL)
  
  if ("BIC-R" %in% Criteria) {
  relaxed = sel.BIC.relaxed(ALL_COEFS,X1,Y2,Z1,pairs,p,K)
  }
  
  ALL_COEFSIter=NULL
  CVMEANIter=CVSEIter=NULL
  cvmeaniter=cv_1se_iter=cv_1.39se_iter=cv_1.96se_iter=NULL
  
  
  if ("CV-OSL" %in% Criteria) {
  
 
  for (rk in ratio_lam1_lam2) {
    
    rk_mat = sparseMatrix(i = 1:(p*(K+1)), j = 1:(p*(K+1)) , x=c(rep(1/rk,p), rep(1,K*(p))))
    
    WIter=1/abs(ALL_COEFS[,which.min(CVMEAN)]+1e-4)

    
    MatDiagWeightsIterInv = sparseMatrix(i = 1:(p*(K+1)), j = 1:(p*(K+1)) ,x= 1/WIter)
    
    calX.Hier.rat = Xx%*%MatDiagWeightsInv%*%rk_mat%*%MatDiagWeightsIterInv
    
  
    mod0         =  clogitLasso(calX.Hier.rat,Y2,strata = pairs,nbfraction = 50,epsilon = ratio)
  

    ################
      
    mod1        =  CV.clogit1(cbind(Y2,Xx%*%MatDiagWeightsInv%*%rk_mat,Y1),NumFolds,mod0$fraction)
    verifI=mod1$verif
  
    CVMEANIter= c(CVMEANIter, mod1$meancv)
    CVSEIter=c(CVSEIter,mod1$secv)
  

    bbeta		 = t(mod0$beta)
    
    
    ALL_COEFSIter = cbind(ALL_COEFSIter, apply(bbeta,2,function(v){as.vector(Matrix(diag(c(rep(1/rk,p), rep(1,K*(p)))), sparse=T)%*%MatDiagWeightsInv%*%MatDiagWeightsIterInv%*%v)}))
    
    
  }
  
  }
  
  if ("BIC" %in% Criteria) {
   bic=refo(vec = ALL_COEFS[,which.min(BIC)],K)
  }
  if ("BIC-R" %in% Criteria) {
    bicR=refo(vec = relaxed[[1]][[which.min(relaxed[[2]])]]$coef.relax,K)
  }
  if ("CV" %in% Criteria || "CV-OSL" %in% Criteria) {
  cvmean= refo(vec = ALL_COEFS[,which.min(CVMEAN)],K)
  
  for (cc in c(1,1.39,1.96)) {
  
  Cond= CVMEAN[which.min(CVMEAN)] + CVSE[which.min(CVMEAN)]*cc
  POsition=which(CVMEAN<=Cond)
  
  if (length(POsition)>1) {
  Comp=apply(ALL_COEFS[,POsition], 2, function(b){
    length(which(b!=0))
  })
  OneSe=POsition[which.min(Comp)]
  }else{OneSe=POsition}
  
  eval(parse(text = paste0("cv_",cc,"se  = refo(vec =ALL_COEFS[,OneSe] ,K)")))
  }
  }
  if ("CV-OSL" %in% Criteria) {
    cvmeaniter= refo(vec = ALL_COEFSIter[,which.min(CVMEANIter)],K)
    
    
    for (cc in c(1,1.39,1.96)) {
    
    Cond= CVMEANIter[which.min(CVMEANIter)] + CVSEIter[which.min(CVMEANIter)]*cc
    POsition=which(CVMEANIter<=Cond)
    
    if (length(POsition)>1) {
    Comp=apply(ALL_COEFSIter[,POsition], 2, function(b){
      length(which(b!=0))
    })
    OneSeiter=POsition[which.min(Comp)]
    }else{OneSeiter=POsition}
    
    eval(parse(text = paste0("cv_",cc,"se_iter  = refo(vec =ALL_COEFSIter[,OneSeiter] ,K)")))
    }
    
  }
  
  bic=bic;bicR=bicR;  cvmean=cvmean;cv_1se=cv_1se;cv_1.39se=cv_1.39se;cv_1.96se=cv_1.96se;cvmeaniter=cvmeaniter;cv_1se_iter=cv_1se_iter;cv_1.39se_iter=cv_1.39se_iter;cv_1.96se_iter=cv_1.96se_iter;
  ALL_COEFS=ALL_COEFS;BIC=BIC;CVMEAN=CVMEAN;CVSE=CVSE;ALL_COEFSR=relaxed[[1]];BICR=relaxed[[2]];ALL_COEFSIter=ALL_COEFSIter;CVMEANIter=CVMEANIter;CVSEIter=CVSEIter
  
return(list(BIC=bic,BIC_R=bicR,CV=cvmean,CV.1se=cv_1se,CV_OSL=cvmeaniter,CV.1se_OSL=cv_1se_iter))
  }
  }




