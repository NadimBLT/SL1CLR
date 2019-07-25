Hierise = function(Feat, ncs)
{
  # block diagonal bg(X1,X2,....,XK) puis 1{C=1},...,1{C=K} puis X, puis 1
  Fin = cbind(Feat, Calise(Feat,ncs))
  return(Fin)
}

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

mise_forme_finale<-function(vect,nstr)
{
  matrix(vect,nrow=nstr,byrow=T)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}



  compute.beta.bar.DELTA = function(coefs.hier, p, N_Str)
  {
    beta.bar        = coefs.hier[1:p]
    DELTA           = matrix(coefs.hier[-(1:p)],nrow = N_Str,byrow = TRUE)
    return(list(beta.bar, DELTA))
  } 

  


sel.BIC.relaxed = function(ALL_COEFS, X, y, Z, S, p, N_Str)
{		
  
  comp.BIC  = function(coefs.hier) 
  {
    DECOMP = compute.beta.bar.DELTA(coefs.hier,p,N_Str)
    return(comp.BIC.relaxed(DECOMP[[1]],DECOMP[[2]], X, y, Z, S, p, N_Str))
  } 
  TEST_ALL  = apply(ALL_COEFS, 2,comp.BIC)
  BICs      = unlist(lapply(TEST_ALL, FUN=function(x){x$Bic}))
  return(list(
    TEST_ALL,
    BICs
  ))	
}


comp.BIC.relaxed = function(beta.bar,DELTA, X, y, Z, S, p, N_Str)
{
  BIC=NULL
  X1=X
  paire=seq(2,nrow(X1),by = 2)
  impaire=seq(1,nrow(X1),by = 2)
  ncs=table(Z)
  calX      = Calise(X1, ncs)
  Ind.p     = 1:p
  Ind.Str   = 1:N_Str
  n         = sum(ncs)
  
  betabar.relax   = rep(0,p)
  deltas.relax    = rep(0,N_Str*p)
  
  ind.hetero      = vector('list')
  for (jj in 1:p)
  {
    ind.hetero[[jj]] = which(DELTA[,jj]!= 0)
  }
  only.common     = which(as.numeric(lapply(ind.hetero, FUN=function(x) {length(x)==0}))==1)
  # positions where only beta.bar is possibly non null
  only.spec       = which(as.numeric(lapply(ind.hetero, FUN=function(x) {length(x)==N_Str}))==1)
  # positions where delta is non null in all strata: these positions can be removed from the overall Feature matrix
  beta.bar.zeros  = which(beta.bar==0)
  # positions where beta.bar is null: they can also be removed from the overall Feature matrix
  mixed           = Ind.p[!(Ind.p%in% only.common | Ind.p%in%only.spec)]
  
  # other positions: for these ones, the overall Feature matrix and the non-null strata have to be retained
  
  if (length( unique(c(only.spec,beta.bar.zeros)))>0)
  {
    Feat            = as.matrix(X1[,-unique(c(only.spec,beta.bar.zeros))])
  }else {Feat = X1}
  ind.remove.calX = NULL
  if (length(only.common)>0) {ind.remove.calX = as.numeric(sapply(only.common, FUN=function(x){x + (0:(N_Str-1))*p}))}
  
  pos.remo = NULL
  if (length(mixed)>0)
  {			
    for (col in mixed)
    {
      Str.remo  = Ind.Str[!Ind.Str%in%ind.hetero[[col]]]
      pos.remo  = c(pos.remo, col+ (Str.remo-1)*p)
    }
  }		
  ind.remove.calX = c(ind.remove.calX,pos.remo)
  if(is.null(ind.remove.calX)==FALSE){
    ind.remove.calX = sort(ind.remove.calX ) 
  }
  
  
  Length.beta.bar  = ncol(Feat)			
  Length.delta     = N_Str*p - length(ind.remove.calX)
  
  if(is.null(ind.remove.calX)==FALSE){
    Features = cbind(Feat, calX[,-ind.remove.calX])
  }else{
   Features = cbind(Feat, calX) 
  }
  
  
  
  if (ncol(Features) >= 1)
  {
     
        if(ncol(Features)>1){
        
        mod0         =  clogitLasso(as.matrix(Features),y,strata = S,nbfraction = 2,epsilon = 1e-10)
        coefs        = (as.matrix(mod0$beta))
        
        for (l in 1:dim(mod0$beta)[1]) {
          XB1=exp(as.vector(as.matrix(Features)[impaire,]%*%mod0$beta[l,]))
          XB2=exp(as.vector(as.matrix(Features)[paire,]%*%mod0$beta[l,]))
          
          XB1[which(XB1==Inf)]=1e+300
          XB1[which(XB1==-Inf)]=-1e+300
          
          XB2[which(XB2==Inf)]=1e+300
          XB2[which(XB2==-Inf)]=-1e+300
          
          XB1[intersect(which(XB2==0),which(XB1==0))]=1e-300
          
          LX=XB1/(XB1+XB2)
          
          LX[which(LX==0)]=1e-300
          
          LL = sum(log(LX))

          
          LL=as.numeric(LL)
          
          BICC= -2*LL + length(which(mod0$beta[l,]!=0))*log(length(y))
          
          BIC=c(BIC,BICC)
        }
     
    
    mod.complex  = ncol(Features)
    
    coef.relaxed = coefs[2,]
    
    coef.temp   = coef.relaxed 
    
    if (Length.beta.bar >= 1) 
    {
      betabar.relax[Ind.p[!Ind.p%in%unique(c(only.spec,beta.bar.zeros))]] = coef.relaxed[1:Length.beta.bar]
      coef.temp   = coef.temp[-(1:Length.beta.bar)]
    }
    if (Length.delta >= 1)
    {
      TEMPO          = 1:(N_Str*p)
      deltas.relax[TEMPO[!TEMPO%in%ind.remove.calX]] = coef.temp[1:Length.delta]
      coef.temp      = coef.temp[-(1:Length.delta)]
    }
    
    BIC1=BIC[2]
    
    
        }else{
          Features=as.matrix(Features)
          mod0 = clogit(y~Features + strata(S),data.frame(y,Features,S))
          BIC1= -2*mod0$loglik[2]+ length(which(mod0$coefficients!=0))*log(length(y))
          coef.relaxed=mod0$coefficients
          
          
          coef.temp   = coef.relaxed 
          
          if (Length.beta.bar >= 1) 
          {
            betabar.relax[Ind.p[!Ind.p%in%unique(c(only.spec,beta.bar.zeros))]] = coef.relaxed[1:Length.beta.bar]
            coef.temp   = coef.temp[-(1:Length.beta.bar)]
          }
          if (Length.delta >= 1)
          {
            TEMPO          = 1:(N_Str*p)
            deltas.relax[TEMPO[!TEMPO%in%ind.remove.calX]] = coef.temp[1:Length.delta]
            coef.temp      = coef.temp[-(1:Length.delta)]
          }
          
          
          
        }
    
  }else {BIC1 = Inf}  

  
  
  coef.relax = c(betabar.relax, deltas.relax)

  
  return(list(
Bic=BIC1,
    coef.relax = coef.relax
  ))
}







sel.BIC.relaxed_IND = function(ALL_COEFS, X, y, S, p)
{		
  
  comp.BIC_IND  = function(coefs.hier) 
  {
    
    return(comp.BIC.relaxed_IND(coefs.hier, X, y, S, p))
    
  } 
  TEST_ALL  = apply(ALL_COEFS, 2,comp.BIC_IND)
  BICs      = unlist(lapply(TEST_ALL, FUN=function(x){x$Bic}))
  return(list(
    TEST_ALL,
    BICs
  ))	
}

comp.BIC.relaxed_IND <- function(coefs.hier, X, y, S, p){
  
  BIC=NULL
  paire=seq(2,nrow(X),by = 2)
  impaire=seq(1,nrow(X),by = 2)
  coef.relaxed=rep(0,p)
  ind=which(coefs.hier!=0)
  if(length(ind)>=1){
    
    if(length(ind)>1){
      X=X[,ind]
      mod0         =  clogitLasso(X,y,strata = rep(1:(length(y)/2),each=2),nbfraction = 2,epsilon = 1e-10)
      for (l in 1:dim(mod0$beta)[1]) {
        XB1=exp(as.vector(X[impaire,]%*%mod0$beta[l,]))
        XB2=exp(as.vector(X[paire,]%*%mod0$beta[l,]))
        
        XB1[which(XB1==Inf)]=1e+300
        XB1[which(XB1==-Inf)]=-1e+300
        
        XB2[which(XB2==Inf)]=1e+300
        XB2[which(XB2==-Inf)]=-1e+300
        
        XB1[intersect(which(XB2==0),which(XB1==0))]=1e-300
        
        LX=XB1/(XB1+XB2)
        
        LX[which(LX==0)]=1e-300
        
        LL = sum(log(LX))
        
        LL=as.numeric(LL)
        
        BICC= -2*LL + length(which(mod0$beta[l,]!=0))*log(length(y))
        
        BIC=c(BIC,BICC)
      }
      
      BIC=BIC[2]
      
      coef.relaxed[ind] = mod0$beta[2,]
    }else{
      X=X[,ind]
      X=as.matrix(X)
      s=rep(1:(length(y)/2),each=2)
      mod0 = clogit(y~X + strata(s),data.frame(y,X,s))
      BIC= -2*mod0$loglik[2]+ length(which(mod0$coefficients!=0))*log(length(y))
      coef.relaxed[ind]=mod0$coefficients
    }
    
  }else{BIC=Inf}
  
  return(list(
    Bic=BIC,
    coef.relax = coef.relaxed
  ))
  
}



