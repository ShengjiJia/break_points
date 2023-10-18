library(glmnet)
library(splines)
library(mgcv)

################Real Data 2
data <- read.csv("C:/Users/Acer/Desktop/My Document/Documents/Research/Projects/change points(series)/Australian CPI data.csv",stringsAsFactors = FALSE)
View(data)
x=data$Food
x=ts(x,start=c(1981,1),end=c(2023,2),frequency=4)
y=data$All.groups.CPI
y=ts(y,start=c(1981,1),end=c(2023,2),frequency=4)
n=length(x)
t=(1:n)/n
X=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1)) X[i,(i+1):n]=0
X1=X[,-1]
grid1=10^seq(1,-1,length=50)    #grid of lambda1
grid2=10^seq(1,-1,length=50)   #grid of lambda2
X2=X1
for(i in 1:n) X2[i,]=X1[i,]*x[i]
X3=cbind(X1,X2)
CV=rep(0,50)          #CV score for each fixed lambda1
bestlam=rep(0,50)    #best lambda2 for each fixed lambda1
for(k in 1:50){
  yy=gam(y~s(t)+s(t,by=x),sp=c(grid1[k],grid1[k]))$residuals
  xx=NULL
  for(i in 1:(2*n-2)){
    xx=cbind(xx,gam(X3[,i]~s(t)+s(t,by=x),sp=c(grid1[k],grid1[k]))$residuals)
  }
  cv.out=cv.glmnet(xx,yy,alpha=1,intercept=FALSE,lambda=grid2,standardize=FALSE)
  CV[k]=min(cv.out$cvm)
  bestlam[k]=cv.out$lambda.min    #best lambda2
}
lambda1=grid1[which.min(CV)]     #best lambda1
lambda2=bestlam[which.min(CV)]     #best lambda2
yy=gam(y~s(t)+s(t,by=x),sp=c(lambda1,lambda1))$residuals
xx=NULL
for(i in 1:(2*n-2)){
  xx=cbind(xx,gam(X3[,i]~s(t)+s(t,by=x),sp=c(lambda1,lambda1))$residuals)
}
beta=glmnet(xx,yy,alpha=1,intercept=FALSE,lambda=lambda2,standardize=FALSE)$beta
model=gam(y-as.vector(X3%*%beta)~s(t)+s(t,by=x),sp=c(lambda1,lambda1))
fit=model$fitted.values+as.vector(X3%*%beta)    #fitted curve

par(mfrow=c(3,2))
plot(as.vector(x)~t, type="l", ylab="x",main="food index")
plot(as.vector(y)~t, type="l", ylab="y",main="all-group index")
lines(fit~t, lty=2, lwd=2, col="red")
fit1=predict(model,newdata=data.frame(t=(1:n)/n,x=rep(0,n)))+as.vector(X1%*%beta[1:(n-1)])
plot(fit1~t,type="l",ylab="f1",main="f1 (proposed method)")
fit2=predict(model,newdata=data.frame(t=(1:n)/n,x=rep(1,n)))+as.vector(cbind(X1,X1)%*%beta)-fit1
plot(fit2~t,type="l",ylab="f2",main="f2 (proposed method)")
#points((which(beta!=0)[-which(diff(which(beta!=0))<=3)]-(n-1))/n,fit2[which(beta!=0)[-which(diff(which(beta!=0))<=3)]-(n-1)], pch=2)
plot(as.vector(y)-fit~t,type="l", ylab="residual", main="Pearson's residuals")
acf(as.vector(y)-fit, main="autocorrelation coefficients")
