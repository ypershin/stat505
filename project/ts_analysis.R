library(tseries)
library(xts)
library(forecast)
library(plyr)

setwd("C:/Users/me/Google Drive/Learning/STAT505/project")


fn_calc <- function(x, lags_=5) {
  mu <- mean(x)
  denom = sum(sapply(x, function(x) { (x-mu)^2 }))
  
  
  rho <- c(); phi <- c()
  for(k in 1:lags_) {
    nom = 0
    for(i in 1:(length(x)-k)) {
      nom = nom + (x[i]-mu)*(x[i+k]-mu)
    }
    rho[k] =  nom/denom
    if(k==1) 
      phi[11] <- rho[1]
    else {
      sum_ = c(0,0)
      n=k-1
      for(j in 1:n) {
        sum_[1] = sum_[1] + phi[n*10+j]*rho[k-j]
        sum_[2] = sum_[2] + phi[n*10+j]*rho[j]
      }
      phi[k*11] = (rho[k]-sum_[1])/(1-sum_[2])
      
      for(j in 1:n) {
        phi[k*10+j] = phi[n*10+j] - phi[k*11]*phi[n*10+(k-j)]
      }
    }
    
  }
  phi_ <- c()
  for(k in 1:length(rho)) {
    # print(sprintf("rho[%d] %0.3f, phi[%d%d] %0.3f",k,rho[k],k,k,phi[k*11]))
    phi_[k] <- phi[k*11]
  }
  return(list(rho=rho,phi=phi_))
}

# plot
fn_plot <- function(x,ttl="") {
  
  # lag 1 difference of z
  xd <- diff(x)
  cf <- fn_calc(x)
  cfd <- fn_calc(xd)
  
  # Augmented Dickey-Fuller test
  a <- adf.test(x); ad <- adf.test(xd);
  
  par(mfrow=c(3,2))
  
  plot.ts(x,main=sprintf("%s, adf %0.2f",ttl,a$p.value),xlab=sprintf("%d observations",length(x)),ylab='')
  # plot.ts(x,main=sprintf("%s",ttl),xlab=sprintf("%d observations",length(x)),ylab='')
  abline(h=mean(x),col="green")
  lines(fitted(tslm(x~trend)),col="red")
  
  plot.ts(xd,main=sprintf("diff(lag=1), adf %0.2f",ad$p.value),xlab=" ",ylab='')
  # plot.ts(xd,main=sprintf("diff(lag=1)"),xlab=" ",ylab='')
  abline(h=mean(xd),col="green")
  lines(fitted(tslm(xd~trend)),col="red")
  
  acf(x,main=sprintf("rho[1..5] %0.2f,%0.2f,%0.2f,%0.2f,%0.2f",cf$rho[1],cf$rho[2],cf$rho[3],cf$rho[4],cf$rho[5]))
  acf(xd,main=sprintf("rho[1..5] %0.2f,%0.2f,%0.2f,%0.2f,%0.2f",cfd$rho[1],cfd$rho[2],cfd$rho[3],cfd$rho[4],cfd$rho[5]))
  # acf(x)
  # acf(xd)
  
  pacf(x,main=sprintf("phi[1..5] %0.2f,%0.2f,%0.2f,%0.2f,%0.2f",cf$phi[1],cf$phi[2],cf$phi[3],cf$phi[4],cf$phi[5]))
  pacf(xd,main=sprintf("phi[1..5] %0.2f,%0.2f,%0.2f,%0.2f,%0.2f",cfd$phi[1],cfd$phi[2],cfd$phi[3],cfd$phi[4],cfd$phi[5]))
  # pacf(x)
  # pacf(xd)
  
  
  par(mfrow=c(1,1))
  
  
}





file_name = "CAD_USD"

df <- read.csv(sprintf("%s.csv",file_name),stringsAsFactors = F)
# head(df)

# x <- ts(df$m1atmiv,start=as.Date(df$date[1],"%Y-%m-%d"),frequency=1)
x <- ts(df$USD,start=as.Date(sprintf("%s-01",df$date[1]),"%Y-%m-%d"),frequency=1)
# head(x,10)


fn_plot(x,file_name)
fn_plot(diff(x))


MAX_ORDER = 5

aic <- NULL
for(p in 0:MAX_ORDER) {
  for(d in 1:1) {
    for(q in 0:MAX_ORDER) {
      tryCatch ({
        x.fit <- arima(x,c(p,d,q))
        aic <- rbind(aic,data.frame(p=p,d=d,q=q,value=x.fit$aic))
      }, error = function(e) {
        print(paste(p,d,q,"non-stationary"))
      })
    }
  }
}

aa <- arrange(aic,value)
head(aa)


model <- as.numeric(aa[1,-4])

x.fit <- arima(x,model)
x.fit
tsdiag(x.fit)

par(mfrow=c(2,1))
pacf(x.fit$residuals)

qqnorm(x.fit$residuals)
par(mfrow=c(1,1))

# cf <- as.vector(x.fit$coef)
# z <- arima.sim(n=length(x),list(order=model,ar=cf[1:2],ma=cf[3:6]))
# fn_plot(z,paste(model,collapse='-'))





