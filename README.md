# SL1CLR
## Stratified L1-penalized conditional logistic regression models for matched case control studies

Methods CondLogist_DataSharedLasso, CondLogist_RefLasso, CondLogist_IndepLasso and CondLogist_PooledLasso, developped by Ballout, Garcia and Viallon [1], are used to estimate several L1-penalized conditional logistic regression models defined on several strata for matched case control studies. Strata can be defined according to specific covariates, but also disease subtypes. Functions corresponding to these methods are available in CondLogist_DataSharedLasso.R, CondLogist_RefLasso.R, CondLogist_IndepLasso.R and CondLogist_PooledLasso.R respectively).

**CondLogist_PooledLasso** estimates the models in a pooled way by fitting L1-penalized conditional logistic regression after pooling all strata together.

**CondLogist_IndepLasso** estimates the models in an independent way by fitting L1-penalized conditional logistic regression on each stratum separately.

**CondLogist_RefLasso** estimates the models in a semi-joint way by fitting L1-penalized conditional logistic regression after choosing a reference stratum a priori and using an additive decomposition of parameters.

**CondLogist_DataSharedLasso** estimates the models in a joint way by fitting L1-penalized conditional logistic regression after using the data shared lasso method (Ollier and Viallon 2017) [4].

See Illustr.R for a simple example illustrating the use of these functions.

See (https://arxiv.org/abs/1901.01583) [1] for more details.
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
* **matchSet**   : Categorical variable indicating sets of matched case-control individuals (only pairs are supported).
* **method**        : Character string, specifies the tuning parameter selection method to be used. Choices are "BIC", "BIC-R", "CV", and/or "CV-OSL".  
"BIC" :  specifies the **B**ayesian **I**nformation **C**riterion;  
"BIC-R":  specifies the **B**ayesian **I**nformation **C**riterion, adapting the Hybrid-OLS idea (Efron and others 2004) [3];  
"CV"  :  specifies the **C**ross **V**alidation technique;  
"CV-OSL"  :  specifies the **C**ross **V**alidation technique following the ideas of    **O**ne **S**tep **L**asso (Buhlmann and Meier 2008) [2]
* **ref**      : This is for the "CondLogist_RefLasso" method, and allows the user to choose the reference stratum.

#### Value
* **BIC**         : Matrix of parameters obtained for the lambda value that minimized the BIC.    
* **BIC-R**       : Matrix of parameters obtained for the lambda value that minimized the BIC-R.   
* **CV**       : Matrix of parameters obtained for the lambda value that minimized the mean cross validation error.   
* **CV.1se**       : Matrix of parameters obtained for the lambda value that have the error which is within 1 standard error of the minimum of mean cross validation error.  
* **CV-OSL**       : Matrix of parameters obtained for the lambda value that minimized the mean cross validation error in **O**ne **S**tep **L**asso.   
* **CV.1se-OSL**       : Matrix of parameters obtained for the lambda value that have the error which is within 1 standard error of the minimum of mean cross validation error in **O**ne **S**tep **L**asso.

**ALL outputs are matrices, of dimension p x K, where p is the number of variables and K is the number of strata.**  



## References

[1] Ballout, Nadim, Garcia, Cedric, et Viallon, Vivian. Sparse estimation for case-control studies with multiple subtypes of cases. arXiv preprint arXiv:1901.01583, 2019.

[2] Bühlmann, Peter; Meier, Lukas. Discussion: One-step sparse estimates in nonconcave penalized likelihood models. Ann. Statist. 36 (2008), no. 4, 1534--1541.

[3] Efron, Bradley; Hastie, Trevor; Johnstone, Iain; Tibshirani, Robert. Least angle regression. Ann. Statist. 32 (2004), no. 2, 407--499.

[4] Ollier, Edouard; Viallon, Vivian. Regression modelling on stratified data with the lasso. Biometrika 104 (2017), mo. 1, 83–96.
