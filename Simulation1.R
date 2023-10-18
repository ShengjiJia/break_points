library(glmnet)
library(splines)

################Simulation 1
set.seed(1823)
n=200
x=(1:n)/n
X=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1)) X[i,(i+1):n]=0
H=X[,-1]               #matrix H
###########Case 1: discontinuous exp function
signal=exp(x)+1*(x>=0.25)-1*(x>=0.75)
error1=rep(0,200)           #MSE of proposed method
Error1=rep(0,200)           #MSE of smoothing splines
tau1=rep(0,200)              #correct detection of jump point tau1=0.25
Tau1=rep(0,200)             #tau2=0.75
smooth1=rep(0,200)       #automatically select a smooth function
for(j in 1:200){
  CV=rep(0,9)          #CV score for each fixed lambda1(i.e. df)
  grid=10^seq(1,-3,length=100)   #grid of lambda2
  for(k in 1:9){
    y=signal+rnorm(n,mean=0,sd=0.25)
    yy=y-smooth.spline(x,y,df=k+1)$y       #df=2:10
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
  lambda2=grid[min(which(model$cvm<model$cvup[which.min(model$cvm)]))]    #best lambda2
  beta=glmnet(xx,yy,alpha=1,intercept=FALSE,lambda=lambda2,standardize=FALSE)$beta
  fit1=as.vector(H%*%beta)      #fitted piecewise-constant
  fit=smooth.spline(x,y-fit1,df=bestdf)$y+fit1  
  error1[j]=sum((signal-fit)^2)/n
  Error1[j]=sum((signal-smooth.spline(x,y,cv=TRUE,all.knots=TRUE)$y)^2)/n
  if(length(which(beta!=0))){
    tau1[j]=(min(which(beta!=0)-0.25*n)<=3)
    Tau1[j]=(min(which(beta!=0)-0.75*n)<=3)
  }
  smooth1[j]=(sum(abs(beta))==0)
  if(j==1){
    par(mfrow=c(2,2))
    plot(y~x,pch=20, col=8,main="case I")
    lines(x,signal,col="red")
    lines(x,fit,lwd=2,lty=2)
    lines(x,smooth.spline(x,y,cv=TRUE,all.knots=TRUE)$y,col="blue")
  }
}
###########Case 2: discontinuous cos function
signal=cos(2*pi*x)+1*(x>=0.25)-1*(x>=0.75)
error2=rep(0,200)
Error2=rep(0,200)
tau2=rep(0,200)             
Tau2=rep(0,200)             
smooth2=rep(0,200)       
for(j in 1:200){
  CV=rep(0,9)          #CV score for each fixed lambda1(i.e. df)
  grid=10^seq(1,-3,length=100)   #grid of lambda2
  for(k in 1:9){
    y=signal+rnorm(n,mean=0,sd=0.25)
    yy=y-smooth.spline(x,y,df=k+1)$y       #df=2:10
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
  lambda2=grid[min(which(model$cvm<model$cvup[which.min(model$cvm)]))]    #best lambda2
  beta=glmnet(xx,yy,alpha=1,intercept=FALSE,lambda=lambda2,standardize=FALSE)$beta
  fit1=as.vector(H%*%beta)      #fitted piecewise-constant
  fit=smooth.spline(x,y-fit1,df=bestdf)$y+fit1  
  error2[j]=sum((signal-fit)^2)/n
  Error2[j]=sum((signal-smooth.spline(x,y,cv=TRUE,all.knots=TRUE)$y)^2)/n
  if(length(which(beta!=0))){
    tau2[j]=(min(which(beta!=0)-0.25*n)<=3)
    Tau2[j]=(min(which(beta!=0)-0.75*n)<=3)
  }
  smooth2[j]=(sum(abs(beta))==0)
  if(j==1){
    plot(y~x,pch=20, col=8, main="case II")
    lines(x,signal,col="red")
    lines(x,fit,lwd=2,lty=2)
    lines(x,smooth.spline(x,y,cv=TRUE,all.knots=TRUE)$y,col="blue")
  }
}
###########Case 3: smooth exp function
signal=exp(x)
error3=rep(0,200)
Error3=rep(0,200)
tau3=rep(0,200)             
Tau3=rep(0,200)             
smooth3=rep(0,200)       
for(j in 1:200){
  CV=rep(0,9)          #CV score for each fixed lambda1(i.e. df)
  grid=10^seq(1,-3,length=100)   #grid of lambda2
  for(k in 1:9){
    y=signal+rnorm(n,mean=0,sd=0.25)
    yy=y-smooth.spline(x,y,df=k+1)$y       #df=2:10
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
  lambda2=grid[min(which(model$cvm<model$cvup[which.min(model$cvm)]))]    #best lambda2
  beta=glmnet(xx,yy,alpha=1,intercept=FALSE,lambda=lambda2,standardize=FALSE)$beta
  fit1=as.vector(H%*%beta)      #fitted piecewise-constant
  fit=smooth.spline(x,y-fit1,df=bestdf)$y+fit1  
  error3[j]=sum((signal-fit)^2)/n
  Error3[j]=sum((signal-smooth.spline(x,y,cv=TRUE,all.knots=TRUE)$y)^2)/n
  if(length(which(beta!=0))){
    tau3[j]=(min(which(beta!=0)-0.25*n)<=3)
    Tau3[j]=(min(which(beta!=0)-0.75*n)<=3)
  }
  smooth3[j]=(sum(abs(beta))==0)
  if(j==1){
    plot(y~x,pch=20, col=8, main="case III")
    lines(x,signal,col="red")
    lines(x,fit,lwd=2,lty=2)
    lines(x,smooth.spline(x,y,cv=TRUE,all.knots=TRUE)$y,col="blue")
  }
}
###########Case 4: smooth cos function
signal=cos(2*pi*x)
error4=rep(0,200)
Error4=rep(0,200)
tau4=rep(0,200)             
Tau4=rep(0,200)             
smooth4=rep(0,200)
for(j in 1:200){
  CV=rep(0,9)          #CV score for each fixed lambda1(i.e. df)
  grid=10^seq(1,-3,length=100)   #grid of lambda2
  for(k in 1:9){
    y=signal+rnorm(n,mean=0,sd=0.25)
    yy=y-smooth.spline(x,y,df=k+1)$y       #df=2:10
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
  lambda2=grid[min(which(model$cvm<model$cvup[which.min(model$cvm)]))]    #best lambda2
  beta=glmnet(xx,yy,alpha=1,intercept=FALSE,lambda=lambda2,standardize=FALSE)$beta
  fit1=as.vector(H%*%beta)      #fitted piecewise-constant
  fit=smooth.spline(x,y-fit1,df=bestdf)$y+fit1  
  error4[j]=sum((signal-fit)^2)/n
  Error4[j]=sum((signal-smooth.spline(x,y,cv=TRUE,all.knots=TRUE)$y)^2)/n
  if(length(which(beta!=0))){
    tau4[j]=(min(which(beta!=0)-0.25*n)<=3)
    Tau4[j]=(min(which(beta!=0)-0.75*n)<=3)
  }
  smooth4[j]=(sum(abs(beta))==0)
  if(j==1){
    plot(y~x,pch=20, col=8, main="case IV")
    lines(x,signal,col="red")
    lines(x,fit,lwd=2,lty=2)
    lines(x,smooth.spline(x,y,cv=TRUE,all.knots=TRUE)$y,col="blue")
  }
}
##########################table 1
cbind(c(mean(error1),mean(Error1),mean(error2),mean(Error2),mean(error3),mean(Error3),mean(error4),mean(Error4)),
      c(mean(smooth1),1,mean(smooth2),1,mean(smooth3),1,mean(smooth4),1),
      c(mean(tau1),0,mean(tau2),0,mean(tau3),0,mean(tau4),0),
      c(mean(Tau1),0,mean(Tau2),0,mean(Tau3),0,mean(Tau4),0))