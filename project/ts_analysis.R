library(tseries)
library(xts)
library(forecast)

setwd("C:/Users/Ypershin/Google Drive/Learning/STAT505/project")

k=1

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
  # a <- adf.test(x); ad <- adf.test(xd);
  
  par(mfrow=c(3,2))
  
  # plot.ts(x,main=sprintf("%s, adf %0.2f",ttl,a$p.value),xlab=sprintf("%d observations",length(x)),ylab='')
  plot.ts(x,main=sprintf("%s",ttl),xlab=sprintf("%d observations",length(x)),ylab='')
  abline(h=mean(x),col="green")
  lines(fitted(tslm(x~trend)),col="red")
  
  # plot.ts(xd,main=sprintf("diff(lag=1), adf %0.2f",ad$p.value),xlab=" ",ylab='')
  plot.ts(xd,main=sprintf("diff(lag=1)"),xlab=" ",ylab='')
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





file_name = "OPT-FB.csv"

df <- read.csv(file_name,stringsAsFactors = F)
# head(df)

x <- ts(df$m1atmiv,start=as.Date(df$date[1],"%Y-%m-%d"),frequency=1)
# head(x,10)


fn_plot(x,"FB_IV")


L_MAX = 10

aic <- NULL
for(p in 0:L_MAX) {
  for(d in 0:1) {
    for(q in 0:L_MAX) {
      x.fit <- arima(x,c(p,d,q))
      aic <- rbind(aic,data.frame(p=p,d=d,q=q,value=x.fit$aic))
    }
  }
}

aa <- arrange(aic,value)
print(aa[1:5,])

nn=1

model <- c(aa$p[nn],aa$d[nn],aa$q[nn])

x.fit <- arima(x,model)
tsdiag(x.fit)
pacf(x.fit$residuals)

qqnorm(x.fit$residuals)



# z <- arima.sim(n=length(x),list(order=model,ma=c(-.14,-.07,-.11,0,.04,-.01,-.05,-.10)))
# fn_plot(z,"ARIMA(0,1,8)")



