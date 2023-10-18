library(glmnet)
library(splines)
library(mgcv)

################Simulation 3
set.seed(23)
n=200
t=(1:n)/n
X=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1)) X[i,(i+1):n]=0
X1=X[,-1]
grid1=10^seq(1,-2,length=10)     #grid of lambda1
grid2=10^seq(1,-3,length=100)   #grid of lambda2
###########Scenario 1: discontinuous coefficient f1
f1=exp(t)+1*(t>=0.5)
f2=cos(2*pi*t)
error1=rep(0,100)            #prediction error for proposed method
Error1=rep(0,100)           #estimation error of f1 for proposed method
error2=rep(0,100)            #prediction error for smoothing spline
Error2=rep(0,100)           #estimation error of f1 for smoothing spline
tau=rep(0,100)                  #correct detection of jump point tau=0.5
smooth=rep(0,100)           #automatically select a smooth function
for(j in 1:100){
  x=rnorm(n=200,mean=0,sd=1)
  signal=f1+f2*x
  y=signal+rnorm(n=200,mean=0,sd=0.25)
  X2=X1
  for(i in 1:n) X2[i,]=X1[i,]*x[i]
  X3=cbind(X1,X2)
  CV=rep(0,10)          #CV score for each fixed lambda1
  bestlam=rep(0,10)    #best lambda2 for each fixed lambda1
  for(k in 1:10){
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
  error1[j]=sum((signal-fit)^2)/n
  Error1[j]=sum((predict(model,newdata=data.frame(t=(1:n)/n,x=rep(0,n)))+as.vector(X1%*%beta[1:(n-1)])-f1)^2)/n
  if(length(which(beta[1:(n-1)]!=0))) tau[j]=(min(which(beta[1:(n-1)]!=0)-0.5*n)<=3)
  smooth[j]=(sum(abs(beta[1:(n-1)]))==0)
  Model=gam(y~s(t)+s(t,by=x))    #standard varying coefficient model
  error2[j]=sum((signal-Model$fitted.values)^2)/n
  Error2[j]=sum((predict(Model,newdata=data.frame(t=(1:n)/n,x=rep(0,n)))-f1)^2)/n
  if(j==1){
    par(mfrow=c(3,2))
    plot(predict(model,newdata=data.frame(t=(1:n)/n,x=rep(0,n)))+as.vector(X1%*%beta[1:(n-1)])~t,type="l",col="blue",lty=2,lwd=2,ylab="f1",main="f1 (proposed method)")
    lines(f1~t,col="red")
    plot(predict(Model,newdata=data.frame(t=(1:n)/n,x=rep(0,n)))~t,type="l",col="blue",lty=2,lwd=2,ylab="f1",main="f1 (smoothing spline)")
    lines(f1~t,col="red")
    plot(model,se=FALSE,rug=FALSE,col="blue", lty=2, lwd=2, ylab="f2",select=2,main="f2 (proposed method)")
    lines(f2~t,col="red")
    plot(Model,se=FALSE,rug=FALSE,col="blue", lty=2, lwd=2, ylab="f2",select=2,main="f2 (smoothing spline)")
    lines(f2~t,col="red")
    plot(y~t, xlim=c(0.4,0.6), ylab="y",main="proposed method")
    lines(t,f1+f2*x,col="red")
    lines(t,fit,lty=2,lwd=2,col="blue")
    plot(y~t, xlim=c(0.4,0.6), ylab="y",main="smoothing spline")
    lines(t,f1+f2*x,col="red")
    lines(t,Model$fitted.values,lty=2,lwd=2,col="blue")
  }
}

###########Scenario 2: continuous coefficient f1
f1=exp(t)
f2=cos(2*pi*t)
error3=rep(0,100)            #prediction error for proposed method
Error3=rep(0,100)           #estimation error of f1 for proposed method
error4=rep(0,100)            #prediction error for smoothing spline
Error4=rep(0,100)           #estimation error of f1 for smoothing spline
Tau=rep(0,100)                  #correct detection of jump point tau=0.5
Smooth=rep(0,100)           #automatically select a smooth function
for(j in 1:100){
  x=rnorm(n=200,mean=0,sd=1)
  signal=f1+f2*x
  y=signal+rnorm(n=200,mean=0,sd=0.25)
  X2=X1
  for(i in 1:n) X2[i,]=X1[i,]*x[i]
  X3=cbind(X1,X2)
  CV=rep(0,10)          #CV score for each fixed lambda1
  for(k in 1:10){
    yy=gam(y~s(t)+s(t,by=x),sp=c(grid1[k],grid1[k]))$residuals
    xx=NULL
    for(i in 1:(2*n-2)){
      xx=cbind(xx,gam(X3[,i]~s(t)+s(t,by=x),sp=c(grid1[k],grid1[k]))$residuals)
    }
    cv.out=cv.glmnet(xx,yy,alpha=1,intercept=FALSE,lambda=grid2,standardize=FALSE)
    CV[k]=min(cv.out$cvm)
  }
  lambda1=grid1[which.min(CV)]     #best lambda1
  yy=gam(y~s(t)+s(t,by=x),sp=c(lambda1,lambda1))$residuals
  xx=NULL
  for(i in 1:(2*n-2)){
    xx=cbind(xx,gam(X3[,i]~s(t)+s(t,by=x),sp=c(lambda1,lambda1))$residuals)
  }
  cv.out=cv.glmnet(xx,yy,alpha=1,intercept=FALSE,lambda=grid2,standardize=FALSE)
  lambda2=grid2[min(which(cv.out$cvm<cv.out$cvup[which.min(cv.out$cvm)]))]    #best lambda2
  beta=glmnet(xx,yy,alpha=1,intercept=FALSE,lambda=lambda2,standardize=FALSE)$beta
  model=gam(y-as.vector(X3%*%beta)~s(t)+s(t,by=x),sp=c(lambda1,lambda1))
  fit=model$fitted.values+as.vector(X3%*%beta)    #fitted curve
  error3[j]=sum((signal-fit)^2)/n
  Error3[j]=sum((predict(model,newdata=data.frame(t=(1:n)/n,x=rep(0,n)))+as.vector(X1%*%beta[1:(n-1)])-f1)^2)/n
  if(length(which(beta[1:(n-1)]!=0))) Tau[j]=(min(which(beta[1:(n-1)]!=0)-0.5*n)<=3)
  Smooth[j]=(sum(abs(beta[1:(n-1)]))==0)
  Model=gam(y~s(t)+s(t,by=x))    #standard varying coefficient model
  error4[j]=sum((signal-Model$fitted.values)^2)/n
  Error4[j]=sum((predict(Model,newdata=data.frame(t=(1:n)/n,x=rep(0,n)))-f1)^2)/n
  if(j==1){
    par(mfrow=c(3,2))
    plot(predict(model,newdata=data.frame(t=(1:n)/n,x=rep(0,n)))+as.vector(X1%*%beta[1:(n-1)])~t,type="l",col="blue",lty=2,lwd=2,ylab="f1",main="f1 (proposed method)")
    lines(f1~t,col="red")
    plot(predict(Model,newdata=data.frame(t=(1:n)/n,x=rep(0,n)))~t,type="l",col="blue",lty=2,lwd=2,ylab="f1",main="f1 (smoothing spline)")
    lines(f1~t,col="red")
    plot(model,se=FALSE,rug=FALSE,col="blue", lty=2, lwd=2, ylab="f2",select=2,main="f2 (proposed method)")
    lines(f2~t,col="red")
    plot(Model,se=FALSE,rug=FALSE,col="blue", lty=2, lwd=2, ylab="f2",select=2,main="f2 (smoothing spline)")
    lines(f2~t,col="red")
    plot(y~t, xlim=c(0.4,0.6), ylab="y",main="proposed method")
    lines(t,f1+f2*x,col="red")
    lines(t,fit,lty=2,lwd=2,col="blue")
    plot(y~t, xlim=c(0.4,0.6), ylab="y",main="smoothing spline")
    lines(t,f1+f2*x,col="red")
    lines(t,Model$fitted.values,lty=2,lwd=2,col="blue")
  }
}

##############show table
cbind(c(mean(error1),mean(error2),mean(error3),mean(error4)),
      c(mean(Error1),mean(Error2),mean(Error3),mean(Error4)),
      c(mean(smooth),1,mean(Smooth),1),
      c(mean(tau),0,mean(Tau),0))
