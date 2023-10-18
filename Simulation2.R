library(glmnet)
library(splines)

################Simulation 2
n=200
X=matrix(1,nrow=n,ncol=n)
for(i in 1:(n-1)) X[i,(i+1):n]=0
x1=X[,-1]
beta=rep(0,n-1)
beta[50]=1
beta[150]=-1
######case 1: null distribution
set.seed(18)
statistic=rep(0,500)
for(k in 1:500){
  x=sort(runif(n,min=0,max=1))
  H=ns(x,knots=x,Boundary.knots = c(0,1))
  #H=cbind(x,cubic(x))
  y=x1%*%beta+rnorm(n,mean=0,sd=0.25)       #null distribution
  grid=10^seq(1,-3,length=100)
  cv.out=cv.glmnet(x1,y,alpha=1,intercept=TRUE,lambda=grid,standardize=FALSE)
  bestlam=cv.out$lambda.min
  fitted=predict(cv.out,newx=x1,s=bestlam)
  residual=y-fitted
  U=0
  for(i in 1:n){
    for(j in 1:n){
      if(j!=i)
        U=U+t(H[i,])%*%H[j,]*residual[i]*residual[j]/n
    }
  }
  R=0
  for(i in 1:n){
    for(j in 1:n){
      if(j!=i)
        R=R+(t(H[i,])%*%H[j,])^2*(residual[i])^2*(residual[j])^2/n/(n-1)
    }
  }
  statistic[k]=abs(U)/sqrt(2*R)           #test statistic
}
par(mfrow=c(2,2))
plot(qnorm(p=0.5+(1:99)/200),quantile(statistic,probs=(1:99)/100), xlab="Theoritical Quantiles", ylab="Sample Quantiles",main="QQ-plot")
lines(qnorm(p=0.5+(1:99)/200),qnorm(p=0.5+(1:99)/200),lwd=2)          #QQ-plot 
######case 1: power
alpha=0.05
power=rep(0,7)        #power function
power[1]=mean(statistic>qnorm(1-alpha/2,mean=0,sd=1))           #power under null distribution
gamma=1:6
for(l in 1:6){
  statistic=rep(0,200)
  for(k in 1:200){
    x=sort(runif(n,min=0,max=1))
    H=ns(x,knots=x,Boundary.knots = c(0,1))
    y=x1%*%beta+rnorm(n,mean=0,sd=0.25)+gamma[l]*exp(x)
    grid=10^seq(1,-3,length=100)
    cv.out=cv.glmnet(x1,y,alpha=1,intercept=TRUE,lambda=grid,standardize=FALSE)
    bestlam=cv.out$lambda.min
    fitted=predict(cv.out,newx=x1,s=bestlam)
    residual=y-fitted
    U=0
    for(i in 1:n){
      for(j in 1:n){
        if(j!=i)
          U=U+t(H[i,])%*%H[j,]*residual[i]*residual[j]/n
      }
    }
    R=0
    for(i in 1:n){
      for(j in 1:n){
        if(j!=i)
          R=R+(t(H[i,])%*%H[j,])^2*(residual[i])^2*(residual[j])^2/n/(n-1)
      }
    }
    statistic[k]=abs(U)/sqrt(2*R)
  }
  power[l+1]=mean(statistic>qnorm(1-alpha/2,mean=0,sd=1))       #power under alternatives
}
plot(power~c(0,gamma),xlab="gamma",ylab="power",ylim=c(0,1),main="Power", type="b")        #power plot

######case 2: null distribution
set.seed(85)
statistic=rep(0,500)
for(k in 1:500){
  x=sort(runif(n,min=0,max=1))
  H=ns(x,knots=x,Boundary.knots = c(0,1))
  #H=cbind(x,cubic(x))
  y=x1%*%beta+rnorm(n,mean=0,sd=0.25)       #null distribution
  grid=10^seq(1,-3,length=100)
  cv.out=cv.glmnet(x1,y,alpha=1,intercept=TRUE,lambda=grid,standardize=FALSE)
  bestlam=cv.out$lambda.min
  fitted=predict(cv.out,newx=x1,s=bestlam)
  residual=y-fitted
  U=0
  for(i in 1:n){
    for(j in 1:n){
      if(j!=i)
        U=U+t(H[i,])%*%H[j,]*residual[i]*residual[j]/n
    }
  }
  R=0
  for(i in 1:n){
    for(j in 1:n){
      if(j!=i)
        R=R+(t(H[i,])%*%H[j,])^2*(residual[i])^2*(residual[j])^2/n/(n-1)
    }
  }
  statistic[k]=abs(U)/sqrt(2*R)           #test statistic
}
plot(qnorm(p=0.5+(1:99)/200),quantile(statistic,probs=(1:99)/100), xlab="Theoritical Quantiles", ylab="Sample Quantiles",main="QQ-plot")
lines(qnorm(p=0.5+(1:99)/200),qnorm(p=0.5+(1:99)/200),lwd=2)          #QQ-plot 
######case 2: power
alpha=0.05
power=rep(0,7)        #power function
power[1]=mean(statistic>qnorm(1-alpha/2,mean=0,sd=1))           #power under null distribution
gamma=(1:6)*0.5
for(l in 1:6){
  statistic=rep(0,200)
  for(k in 1:200){
    x=sort(runif(n,min=0,max=1))
    H=ns(x,knots=x,Boundary.knots = c(0,1))
    y=x1%*%beta+rnorm(n,mean=0,sd=0.25)+gamma[l]*cos(4*pi*x)
    grid=10^seq(1,-3,length=100)
    cv.out=cv.glmnet(x1,y,alpha=1,intercept=TRUE,lambda=grid,standardize=FALSE)
    bestlam=cv.out$lambda.min
    fitted=predict(cv.out,newx=x1,s=bestlam)
    residual=y-fitted
    U=0
    for(i in 1:n){
      for(j in 1:n){
        if(j!=i)
          U=U+t(H[i,])%*%H[j,]*residual[i]*residual[j]/n
      }
    }
    R=0
    for(i in 1:n){
      for(j in 1:n){
        if(j!=i)
          R=R+(t(H[i,])%*%H[j,])^2*(residual[i])^2*(residual[j])^2/n/(n-1)
      }
    }
    statistic[k]=abs(U)/sqrt(2*R)
  }
  power[l+1]=mean(statistic>qnorm(1-alpha/2,mean=0,sd=1))       #power under alternatives
}
plot(power~c(0,gamma),xlab="gamma",ylab="power",ylim=c(0,1),main="Power", type="b")        #power plot
