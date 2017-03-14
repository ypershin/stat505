library(tseries)
library(forecast)


options(warn=-1)

descr <- c(
  'W1: Daily Average Number of Truck Manufacturing Defects.',
  'W2: Wolf Yearly Sunspot Numbers from 1700 to 1983.',
  'W3: Blowfly Data.',
  'W4: Monthly Unemployed Females Between Ages 16 and 19\n in the U.S. from January 1961 to December 1985 (in thousands).',
  'W5: Yearly Accidental Death Rate (per 100,000)\n for the State of Pennsylvania Between 1950 and 1984.',
  'W6: Yearly U.S. Tobacco Production from 1871 to 1984\n (in millions of pounds).',
  'W7: Yearly Number of Lynx Pelts Sold by the Hudson_s Bay Company\n in Canada from 1857 to 1911.',
  'W8: Simulated Seasonal Time Series.',
  'W9: Monthly Employment Figures for Males Between Ages 16 and 19\n in the U.S. from January 1971 to December 1981 (in thousands).',
  'W10: Quarterly U.S. Beer Production from First Quarter of 1975\n to the Fourth Quarter of 1982 (in millions of Barrels).'
)

# (start,frequency)
sf <- list(c(1,1),c(1700,1),c(1,1),c(1961,12),c(1950,1),c(1871,1),c(1857,1),c(1,12),c(1971,12),c(1975,4))

# calculate acf and pcf coeffs
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
fn_plot <- function(x,ttl) {

  # lag 1 difference of z
  xd <- diff(x)
  cf <- fn_calc(x)
  cfd <- fn_calc(xd)
  
  # Augmented Dickey-Fuller test
  a <- adf.test(x); ad <- adf.test(xd);
  
  par(mfrow=c(3,2))
  
  plot.ts(x,main=sprintf("%s, adf %0.2f",ttl,a$p.value),xlab=sprintf("%d observations",length(x)),ylab='')
  abline(h=mean(x),col="green")
  lines(fitted(tslm(x~trend)),col="red")

  plot.ts(xd,main=sprintf("diff(lag=1), adf %0.2f",ad$p.value),xlab=" ",ylab='')
  abline(h=mean(xd),col="green")
  lines(fitted(tslm(xd~trend)),col="red")
  
  acf(x,main=sprintf("rho[1..5] %0.2f,%0.2f,%0.2f,%0.2f,%0.2f",cf$rho[1],cf$rho[2],cf$rho[3],cf$rho[4],cf$rho[5]))
  acf(xd,main=sprintf("rho[1..5] %0.2f,%0.2f,%0.2f,%0.2f,%0.2f",cfd$rho[1],cfd$rho[2],cfd$rho[3],cfd$rho[4],cfd$rho[5]))
  
  pacf(x,main=sprintf("phi[1..5] %0.2f,%0.2f,%0.2f,%0.2f,%0.2f",cf$phi[1],cf$phi[2],cf$phi[3],cf$phi[4],cf$phi[5]))
  pacf(xd,main=sprintf("phi[1..5] %0.2f,%0.2f,%0.2f,%0.2f,%0.2f",cfd$phi[1],cfd$phi[2],cfd$phi[3],cfd$phi[4],cfd$phi[5]))
  
  
  par(mfrow=c(1,1))
  

}




transformed_sqrt <- c(2,3); transformed_log <- c(6,7)

# main
for(i in 1:length(descr)) {
  z <- ts(read.table(sprintf('data/W%d.txt',i)),start=sf[[i]][1],frequency=sf[[i]][2])
  ttl <- substr(descr[i],1,25)
  
  fn_plot(z,ttl)
  if(i %in% c(transformed_sqrt,transformed_log)) {
    
    if(i %in% transformed_sqrt) {
      zt <- sqrt(z); ttl <- paste("[SQRT]",ttl)
    }
    else {
      zt <- log(z); ttl <- paste("[LN]",ttl)
    }
      
    fn_plot(zt,ttl)
  }

  if(sf[[i]][2]>1) {
    plot(decompose(z),xlab=ttl,ylab='')
  }
}

