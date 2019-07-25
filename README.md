# SL1CLR
## Stratified L1-penalized conditional logistic regression models for matched case control studies

Methods CondLogist_DataSharedLasso, CondLogist_RefLasso, CondLogist_IndepLasso and CondLogist_PooledLasso, developped by Ballout, Garcia and Viallon *, are used to estimate several L1-penalized conditional logistic regression models defined on several strata for matched case control studies. Functions corresponding to these methods are available in CondLogist_DataSharedLasso.R, CondLogist_RefLasso.R, CondLogist_IndepLasso.R and CondLogist_PooledLasso.R respectively).

**CondLogist_PooledLasso** estimates the models in a pooled way by fitting L1-penalized conditional logistic regression after pooling all strata together.

**CondLogist_IndepLasso** estimates the models in an independent way by fitting L1-penalized conditional logistic regression on each stratum separately.

**CondLogist_RefLasso** estimates the models in a semi-joint way by fitting L1-penalized conditional logistic regression after choosing a reference stratum a priori and using an additive decomposition of parameters.

**CondLogist_DataSharedLasso** estimates the models in a joint way by fitting L1-penalized conditional logistic regression after using the data shared lasso method (Viallon *).

See Illustr.R for a simple example illustrating the use of these functions.

See https:...... for more details.
## Packages required 


```
cLogitLasso, survival and Matrix.
```



## Usage
##### CondLogist_DataSharedLasso, CondLogist_RefLasso, CondLogist_IndepLasso and CondLogist_PooledLasso
#### Arguments
* **X**        : Input matrix, of dimension n x p, where n is the number of observations and p is the number of variables; each row is an observation vector.  
* **y**        : Binary response variable, with 1 for cases and 0 for controls
* **pairs**    : Vector defining the pairs; each pair composed by a case and his matched control.  
* **strata**   : Categorical variable defining the strata of pairs.  
* **method**        : Character string, specifies the tuning parameter selection method to be used. Choices are "BIC", "BIC-R", "CV", and/or "CV-OSL".  
"BIC" :  specifies the **B**ayesian **I**nformation **C**riterion;  
"BIC-R":  specifies the **B**ayesian **I**nformation **C**riterion **R**elaxed;  
"CV"  :  specifies the **C**ross **V**alidation technique;  
"CV-OSL"  :  specifies the **C**ross **V**alidation technique following the ideas of    **O**ne **S**tep **L**asso
* **ref**      : This is for the "CondLogist_RefLasso" method, and allows the user to choose the reference stratum.

#### Value
* **BIC**         : Estimated models returned by the BIC method, models that have the minmum values of BIC.    
* **BIC-R**       : Estimated models returned by the BIC-R method, models that have the minmum values of BIC in Relaxed version.  .   
* **CV**       : Estimated models returned by the CV method, models that have the minmum values of mean cross validation error.   
* **CV.1se**       : Estimated models returned by the CV method, models that have the error which is within 1 standard error of the minimum of mean cross validation error.  
* **CV-OSL**       : Estimated models returned by the CV method in One Step Lasso, models that have the minmum values of mean cross validation error.   
* **CV.1se-OSL**       : Estimated models returned by the CV method in One Step Lasso, models that have the error which is within 1 standard error of the minimum of mean cross validation error.   
**ALL outputs are matrices, of dimension p x K, where p is the number of variables and K is the number of strata.**  
