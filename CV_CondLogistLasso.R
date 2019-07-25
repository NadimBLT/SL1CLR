
CV.clogit <- function(data,K,lambda){
  
  ncs=table(data[,ncol(data)])[-1]*2
  data=data[,-ncol(data)]
  
  FOLDSS=NULL
  
  for (gg in 1:length(ncs)) {
    
    folds <- rep_len(1:K, ncs[gg]/2) 
    set.seed(100)
    folds <- sample(folds,length(folds),replace = FALSE)
    FOLDSS=c(FOLDSS,folds)
  }
  
  folds <- rep(FOLDSS,each=2)
  
  
  paire_All=seq(2,nrow(data),by = 2)
  impaire_All=seq(1,nrow(data),by = 2)
  
  
  cv.error=matrix(NA,nrow = K,ncol = length(lambda))
  
  for(k in 1:K) {
    
    fold <- which(folds == k)
    data.train <- data[-fold,]
    data.test <- data[fold,]
    
    paire.test=seq(2,nrow(data.test),by = 2)
    impaire.test=seq(1,nrow(data.test),by = 2)
    
    mod0         =  clogitLasso(data.train[,-1],data.train[,1],strata = rep(1:(nrow(data.train)/2),each=2),fraction = lambda)
    
    for (l in 1:dim(mod0$beta)[1]) {
      
      XB1=as.vector(data.test[impaire.test,-1]%*%mod0$beta[l,])
      XB2=as.vector(data.test[paire.test,-1]%*%mod0$beta[l,])
      LL = sum(log((exp(XB1)/(exp(XB1)+exp(XB2)))))
      cv.error[k,l]= as.numeric(LL)
      
    }
    
    
  }
  
  meancv=apply(-cv.error, 2, mean)
  secv=apply(-cv.error, 2, sd)/sqrt(K)
  
  return(list(meancv=meancv,secv=secv,cv.error=cv.error,verif=FOLDSS))
}




CV.clogit1 <- function(data,K,lambda){
  ZD=data[,ncol(data)]
  ncs=table(ZD)[-1]*2
  data=data[,-ncol(data)]
  FOLDSS=NULL
  
  for (gg in 1:length(ncs)) {
    
    folds <- rep_len(1:K, ncs[gg]/2) 
    set.seed(100)
    folds <- sample(folds,length(folds),replace = FALSE)
    FOLDSS=c(FOLDSS,folds)
  }
  
  folds <- rep(FOLDSS,each=2)
  
  paire_All=seq(2,nrow(data),by = 2)
  impaire_All=seq(1,nrow(data),by = 2)
  
  
  cv.error=matrix(NA,nrow = K,ncol = length(lambda))
  
  for(k in 1:K) {
    
    fold <- which(folds == k)
    ZDF <- ZD[-fold]
    data.train <- data[-fold,]
    data.test <- data[fold,]
    
    paire.test=seq(2,nrow(data.test),by = 2)
    impaire.test=seq(1,nrow(data.test),by = 2)
    
    
    mod0F         =  clogitLasso(data.train[,-1],data.train[,1],strata = rep(1:(nrow(data.train)/2),each=2),nbfraction = 50)
    mod1F         =  CV.clogit(data = cbind(data.train,ZDF), lambda = mod0F$fraction,K = 10)
    
    WIter=1/abs(mod0F$beta[which.min(mod1F$meancv),]+1e-4)
    MatDiagWeightsIterInv = sparseMatrix(i = 1:ncol(data.train[,-1]), j = 1:ncol(data.train[,-1]) ,x= 1/WIter)
    mod00         =  clogitLasso(data.train[,-1]%*%MatDiagWeightsIterInv,data.train[,1],strata = rep(1:(nrow(data.train)/2),each=2),fraction = lambda)
    
    
    
    for (l in 1:dim(mod00$beta)[1]) {
      
      XB1=as.vector(data.test[impaire.test,-1]%*%as.vector(mod00$beta[l,]%*%MatDiagWeightsIterInv))
      XB2=as.vector(data.test[paire.test,-1]%*%as.vector(mod00$beta[l,]%*%MatDiagWeightsIterInv))
      LL = sum(log((exp(XB1)/(exp(XB1)+exp(XB2)))))
      cv.error[k,l]= as.numeric(LL)
      
    }
    
    
  }
  
  meancv=apply(-cv.error, 2, mean)
  secv=apply(-cv.error, 2, sd)/sqrt(K)
  
  return(list(meancv=meancv,secv=secv,cv.error=cv.error,verif=FOLDSS))
}
