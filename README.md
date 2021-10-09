# o2plsda: Omics data integration with o2plsda
# o2plsda [![Project Status:](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)  [![](https://img.shields.io/badge/devel%20version-0.0.6-green.svg)](https://github.com/guokai8/o2plsda)  ![Code Size:](https://img.shields.io/github/languages/code-size/guokai8/o2plsda)![](https://img.shields.io/badge/license-GPL--3-blue.svg)[![DOI](https://zenodo.org/badge/413478714.svg)](https://zenodo.org/badge/latestdoi/413478714)


## Description
_o2plsda_ provides functions to do O2PLS-DA analysis for mutiple omics integration. The package could use the group information (MCCV) when we select the best paramaters with cross-validation. The algorithm came from "O2-PLS, a two-block (X±Y) latent variable regression (LVR) method with an integral OSC filter" which published by Johan Trygg and Svante Wold at 2003. The package also add PLS-DA function.

## Installation
```{r,eval=FALSE}
library(devtools)
install_github("guokai8/o2plsda")
``` 
## Examples
```{r}
library(o2plsda)
set.seed(123)
# sample * values
X = matrix(rnorm(5000),50,100)
# sample * values
Y = matrix(rnorm(5000),50,100)
rownames(X) <- paste("S",1:50,sep="")
rownames(Y) <- paste("S",1:50,sep="")
colnames(X) <- paste("Gene",1:100,sep="")
colnames(Y) <- paste("Lipid",1:100,sep="")
X = scale(X, scale=T)
Y = scale(Y, scale=T)
## group factor could be omitted if you don't have any group 
group <- rep(c("Ctrl","Treat"),each = 25)
```
Do cross validation with group information
```{r}
set.seed(123)
## nr_folds : cross validation k-fold (suggest 10)
## ncores : parallel paramaters for large datasets
cv <- o2cv(X,Y,1:5,1:3,1:3,group=group,nr_folds = 10)
#####################################
# The best paramaters are nc =  5 , nx =  3 , ny =  3 
#####################################
# The Qxy is  0.08222935  and the RMSE is:  2.030108 
#####################################
```

Then we can do the O2PLS analysis with nc = 5, nx = 3, ny =3. You can also select the best paramaters by looking at the cross validation results.
```{r}
fit <- o2pls(X,Y,5,3,3)
summary(fit)
######### Summary of the O2PLS results #########
### Call o2pls(X, Y, nc= 5 , nx= 3 , ny= 3 ) ###
### Total variation 
### X: 4900 ; Y: 4900  ###
### Total modeled variation ### X: 0.286 ; Y: 0.304  ###
### Joint, Orthogonal, Noise (proportions) ###
#               X     Y
#Joint      0.176 0.192
#Orthogonal 0.110 0.112
#Noise      0.714 0.696
### Variation in X joint part predicted by Y Joint part: 0.906 
### Variation in Y joint part predicted by X Joint part: 0.908 
### Variation in each Latent Variable (LV) in Joint part: 
#      LV1     LV2     LV3     LV4     LV5
#X 181.764 179.595 191.210 152.174 157.819
#Y 229.308 204.829 175.926 173.382 155.934
### Variation in each Latent Variable (LV) in X Orthogonal part: 
#      LV1     LV2     LV3
#X 227.856 166.718 143.602
### Variation in each Latent Variable (LV) in Y Orthogonal part: 
#      LV1     LV2     LV3
#Y 225.833 166.231 157.976

```

Extract the loadings and scores from the fit results

```{r}
Xl <- loadings(fit,loading="Xjoint")
Xs <- scores(fit,score="Xjoint")
plot(fit,type="score",var="Xjoint", group=group)
plot(fit,type="loading",var="Xjoint", group=group,repel=F,rotation=TRUE)
```

Do the OPLSDA based on the O2PLS results
```{r}
res <- oplsda(fit,group, nc=5)
plot(res,type="score", group=group)
vip <- vip(res)
plot(res,type="vip", group = group, repel = FALSE,order=TRUE)
```

## Note
The package is still under development.

## Citation
If you like this package, please contact me for the citation.

## Contact information

For any questions please contact guokai8@gmail.com or https://github.com/guokai8/o2plsda/issues
