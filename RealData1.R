library(glmnet)
library(splines)

################Real Data 1
load("C:/Users/Acer/Desktop/My Document/Documents/Research/Projects/change points(series)/data4.RData")
y=as.vector(t(data4[4,1:566,1]))      #data20[4st stock, days 1-566, open prices]
n=length(y)
x=(1:n)/n
X=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1)) X[i,(i+1):n]=0
H=X[,-1]               #matrix H
CV=rep(0,49)          #CV score for each fixed lambda1(i.e. df)
grid=10^seq(0,-3,length=100)   #grid of lambda2
set.seed(18)
for(k in 1:49){
  yy=y-smooth.spline(x,y,df=k+1)$y       #df=2:50
  xx=NULL
  for(i in 1:(n-1)){
    xx=cbind(xx,H[,i]-smooth.spline(x,H[,i],df=k+1)$y)
  }
  cv.out=cv.glmnet(xx,yy,alpha=1,intercept=FALSE,lambda=grid,standardize=FALSE)
  CV[k]=min(cv.out$cvm)
}
bestdf=which.min(CV)+1       #best df(or lambda1)
yy=y-smooth.spline(x,y,df=bestdf)$y     #fit with best lambda1
xx=NULL
for(i in 1:(n-1)){
  xx=cbind(xx,H[,i]-smooth.spline(x,H[,i],df=bestdf)$y)
}
model=cv.glmnet(xx,yy,alpha=1,intercept=FALSE,lambda=grid,standardize=FALSE)
lambda2=model$lambda.min
beta=glmnet(xx,yy,alpha=1,intercept=FALSE,lambda=lambda2,standardize=FALSE)$beta
fit1=as.vector(H%*%beta)      #fitted piecewise-constant
fit=smooth.spline(x,y-fit1,df=bestdf)$y+fit1  

####################show Figure 5
par(mfrow=c(2,1))
plot(y~x, type="l", main="daily open prices")
plot(fit~x, type="l", main="fitted curve by proposed method")
abline(v= which(beta!=0)[-which(diff(which(beta!=0))<=3)][-c(2,3)]/n, lty=2)
abline(v= which(beta!=0)[-which(diff(which(beta!=0))<=3)][c(2,3)]/n, lty=2, col="red")
