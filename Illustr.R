path_repository="~/"

library(clogitLasso)
library(survival)
library(Matrix)

source(paste0(path_repository,"CondLogist_DataSharedLasso.R"))
source(paste0(path_repository,"CondLogist_RefLasso.R"))
source(paste0(path_repository,"CondLogist_IndepLasso.R"))
source(paste0(path_repository,"CondLogist_PooledLasso.R"))
source(paste0(path_repository,"CV_CondLogistLasso.R"))
source(paste0(path_repository,"BIC-R_CondLogistLasso.R"))

load(paste0(path_repository,"/DataEx.Rdata"))

# n=1000, p=20, K=6, n_k=(400,200,100,100,100,100)

True_Models=Data$True_Models
X=Data$X
y=Data$y
pairs=Data$pairs
matchSet=Data$strata

True_Models
head(X)
head(y)
head(pairs)
head(matchSet)

'
RES_Pooled       = CondLogist_PooledLasso     (X=X,y=y,pairs=pairs,matchSet=matchSet,method = c("BIC","BIC-R"))
RES_Indep        = CondLogist_IndepLasso      (X=X,y=y,pairs=pairs,matchSet=matchSet,method = c("BIC","BIC-R"))
RES_Ref          = CondLogist_RefLasso        (X=X,y=y,pairs=pairs,matchSet=matchSet,method = c("BIC","BIC-R"),ref=1)
RES_DataShared   = CondLogist_DataSharedLasso (X=X,y=y,pairs=pairs,matchSet=matchSet,method = c("BIC","BIC-R"))
'
load(paste0(path_repository,"ResEx.Rdata"))

# Comparing the estimated models with the true models

#ALL = cbind(True_Models,RES_Pooled$BIC, RES_Indep$BIC,RES_Ref$BIC,RES_DataShared$BIC)
ALL = cbind(True_Models,RES_Pooled$BIC_R, RES_Indep$BIC_R,RES_Ref$BIC_R,RES_DataShared$BIC_R)

ALL[which(ALL>6)]=6
ALL[which(ALL< -6)]=-6
RES_DataFrame = data.frame(DELTA=as.numeric(ALL), METH=rep(c("TRUE", "CondLogist_PooledLasso", "CondLogist_IndepLasso","CondLogist_RefLasso", "CondLogist_DatasharedLasso"), each=20*6), VAR=rep(paste0("Var_",1:20), 5*6), STRAT=  rep(rep(paste0("Stratum_",1:6), each=20), 5))
RES_DataFrame$VAR=factor(RES_DataFrame$VAR,levels = paste0("Var_",1:20))
RES_DataFrame$METH=factor(RES_DataFrame$METH,levels = c("TRUE", "CondLogist_PooledLasso", "CondLogist_IndepLasso","CondLogist_RefLasso", "CondLogist_DatasharedLasso"),labels = c("TRUE","Pooled","Indep","Ref","DataShared"))


library(ggplot2)
ggplot(RES_DataFrame,aes(STRAT,VAR)) + geom_tile(aes(fill = as.numeric(as.character(DELTA))),colour = 'black') + 
  scale_fill_gradient2(high="black",mid="white",low="red",midpoint=0) + 
  facet_wrap(~ METH, ncol=5)+ 
  theme(legend.position="left",legend.title=element_blank(),panel.margin.y = unit(0.8, "lines"),
        strip.background=element_blank(),strip.text=element_text(face="bold",size=1.5),
        axis.ticks = element_line(size = 0.1), axis.ticks.length = unit(0.03, "cm"),
        axis.text=element_text(size=10,face = "bold"), 
        axis.title=element_blank(), axis.text.x=element_text(hjust = -0.4, vjust = 0.4, angle=-90),    
        strip.text.x = element_text(size = 12), legend.text=element_text(size=10,hjust = 1.2,vjust=0.8,face="bold"), legend.key.size = unit(0.5, 'lines'))

