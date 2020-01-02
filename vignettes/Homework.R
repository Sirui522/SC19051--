## ------------------------------------------------------------------------
set.seed(12345)
graph_hist = function(size, sigma){
  
  hist(x = sigma * sqrt(-2*log(runif(size))), prob = T,main="Sample from Rayleigh Distribution",xlab = "value") #histogram of the sample points 
  
  abline(v = sigma,lwd=4,col="blue") #the position of mode
}

graph_hist(10000,0.5)         # Histogram when sigma = 0.5

graph_hist(10000,1)        # Histogram when sigma = 1

graph_hist(10000,5)       # Histogram when sigma = 5




## ------------------------------------------------------------------------
set.seed(54321)

ge_MN = function(n, p1){  # sample of mixture normal
    U = runif(n)
    Y = rep(0,n)
    for(i in 1:n){
      if(U[i] <= p1){
        Y[i] = rnorm(1)
      }
      else{
        Y[i] = rnorm(1) + 3
      }
    }
    
    hist(Y, probability = T, main = "Mixture Normal")
    count = 0
    for(i in 1:n){
      if(Y[i] <= 1.5){
        count = count + 1
      }
    }
    print("The conjecture about p1 is ")
    print(count/n)
}


ge_MN(5000, 0.25)
ge_MN(5000, 0.5)
ge_MN(5000, 0.75)


## ------------------------------------------------------------------------
library(matrixcalc)

ge_Wishart = function(n, Sigma){
  n = as.integer(n)
  d = nrow(Sigma)
  if(n <= d+1){
    return("error in dimension")
  }
  
  UL = chol(Sigma)                  
  # Choleski factorization, where UL is upper triangular
  
  MA = matrix(0, nrow = d, ncol = d)
  
  for(i in 1:d){
    MA[i,i] = sqrt(rchisq(1, n-i+1))
    if(i > 1){
      for(j in 1:(i-1)){
        MA[i,j] = rnorm(1)
      } 
    }
  }
  
  return(t(UL) %*% MA %*% t(MA) %*% UL)
}



## ------------------------------------------------------------------------
B = cbind(c(6,3,0),c(3,6,-6),c(0,-6,11))
ge_Wishart(5,B)

## ------------------------------------------------------------------------
N = 10000 ##size of Monte Carlo sample
est_Mon = pi/3*mean(runif(N,0,pi/3))
true_value = integrate(function(t){sin(t)},lower = 0,upper = pi/3)$value
print(paste("The Monte Carlo estimate is",round(est_Mon,3),"and the true value of the integral is",true_value))

## ------------------------------------------------------------------------
f = function(x){
  return(exp(-x)/(1+x^2))
}  ##define the function f

U1 = runif(N,0,1)
U2 = runif(N,0,1)
est_van = mean(c(f(U1),f(U2)))
est_ant = mean(c(f(U1),f(1-U2)))
var_van = var(c(f(U1),f(U2)))
var_ant = var(c(f(U1),f(1-U2)))
reduction = (var_van - var_ant)/(var_van)
print(paste("The vanilla Monte Carlo estimate is",round(est_van,3)," the estimate with antithetic variables is",round(est_ant,3),"and the true value is",round(integrate(f,lower = 0,upper = 1)$value,3)))
print(paste("The variance of vanilla MC is",round(var_van,4),"and the variance of MC with antithetic variables is",round(var_ant,4),"reducing by",round(reduction,3)))


## ------------------------------------------------------------------------
m <- 10000 #number of replicates
k <- 10 #number of strata
N0 <- 50 #number of times to repeat the estimation 
T = numeric(k)
estimates = matrix(0,N0,2)

g = function(t) { exp(-t - log(1+t^2)) * (t > 0) * (t < 1) }
fg = function(t) { g(t) / (exp(-t) / (1 - exp(-1))) }
inv = function(t) { - log(1 - t * (1 - exp(-1))) } ##inverse transform method 

for (i in 1:N0){
  estimates[i,1] = mean(fg(inv(runif(m))))
  for (j in 1:k){
    T[j] = mean(fg(inv(runif(m/k,(j-1)/k,j/k))))
  estimates[i,2] = mean(T) 
  }
}
v = apply(estimates, 2, var)

print(paste("The estimates of no-stratified and stratified are: ",round(apply(estimates, 2, mean),3)))
print(paste("The variances of no-stratified and stratified are: ",v))
print(paste("This represents a",100*round((v[1]-v[2])/v[1],2),"% reduction in variance."))




## ------------------------------------------------------------------------
set.seed(12345)

m = 5e4                    #MCMC sample size
n = 20                     #sample size

est1 = rep(0,m)            #MC estimate coverage probability of t-interval
est2 = rep(0,m)            #MC estimate coverage probability of chi2-interval

for(i in 1:m){
  
  X = rchisq(n, df = 2)    #sample from chi2 
  
  est1[i] = mean(X)- qt(0.975,n-1)*var(X)/sqrt(n) < 2 & mean(X) + qt(0.975,n-1)*var(X)/sqrt(n) > 2
  # the mean of the chi2 is 2
  
  est2[i] = (n-1)*var(X)/qchisq(0.05, df = n-1) > 4
  # the variance of the chi2 is 4
  
  }

mean(est1)             
mean(est2)           

## ------------------------------------------------------------------------
set.seed(54321)
q0 = c(0.025,0.05,0.95,0.975)
n=500
m = 5e4
skew = rep(0,m)

for(i in 1:m){
  sam = rnorm(n)          
  skew[i] = mean((sam - mean(sam))^3)/(mean((sam - mean(sam))^2))^1.5           
  }
  
est1_q = quantile(skew, probs = q0)           
  
est2_q = qnorm(q0, sd = sqrt(6/n) )         

std1 = sqrt( q0*(1-q0)/n/( dnorm(est1_q, sd = sqrt(6/n)) )^2 )
# the standard error of sample quantile estimates 
  
result = cbind(est2_q, est1_q, std1)
colnames(result) = c("theoretical quantile(Normal)","sample quantile","standard error of sample quantile(2.14)")
  
knitr::kable(result)

## ------------------------------------------------------------------------
sk = function(x) { 
  #sample skewness 
  m1 = mean(x) 
  m3 = mean((x - m1)^3) 
  m2 = mean((x - m1)^2) 
  return( m3 / m2^1.5 ) 
} 

alpha = seq(0,3,0.5)
N = 5000  ##number of iterations
n = 30  ##sample size
sl = 0.05  ##significance level
result1 = rep(0,length(alpha))

for (j in 1:length(alpha)) {           
  cv = qnorm(1-sl/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))  ##critical value for the skewness test
  skvector1 = rep(0,N)
  
  for (i in 1:N) {           
    x = rbeta(n, shape1 = alpha[j], shape2 = alpha[j]) 
    if (abs(sk(x)) >= cv) skvector1[i] = 1
    else skvector1[i] = 0
  } 
  
  result1[j] = mean(skvector1) 
}

plot(alpha, result1, type = "b",xlab = "a", ylab = "power",ylim = c(0,1),main = "Beta(a,a):significance level = 0.05") 
M = rbind(alpha, result1)
rownames(M) = c("alpha","power")
colnames(M) = 1:length(alpha)
knitr::kable(M)

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
sk = function(x) { 
  #sample skewness 
  m1 = mean(x) 
  m3 = mean((x - m1)^3) 
  m2 = mean((x - m1)^2) 
  return( m3 / m2^1.5 ) 
} 

freedom = seq(1,10,1)
N = 5000  ##number of iterations
n = 30  ##sample size
sl = 0.05  ##significance level
result2 = rep(0,length(freedom))

for (j in 1:length(freedom)) {           
  cv = qnorm(1-sl/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))  ##critical value for the skewness test
  skvector2 = rep(0,N)
  
  for (i in 1:N) {           
    x = rt(n, df = freedom[j])  
    if (abs(sk(x)) >= cv) skvector2[i] = 1
    else skvector2[i] = 0
  } 
  
  result2[j] = mean(skvector2) 
}

plot(freedom, result2, type = "b",xlab = "fr", ylab = "power",ylim = c(0,1),main = "t(fr): significance level = 0.05") 
M = rbind(freedom, result2)
rownames(M) = c("degree of fr","power")
colnames(M) = 1:length(freedom)
knitr::kable(M)

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
sl = c(0.01,0.05,0.1,0.2)
N = 3000 ##number of iterations
n = 30 ##sample size
result3 = rep(0,length(sl))

for(j in 1:length(sl)){
    
    t1e_vector = rep(0,n) 

    for(i in 1:N){
      x = rchisq(n, df = 1)                                              
      q = qt(1-sl[j]/2, df = n-1)                                      
      t1e_vector[i] = mean(abs(mean(x)-1)>q*var(x)/sqrt(n))  ## the mean of chi^2(1) is 1
    }
    
    result3[j] = mean(t1e_vector)
  }

result3 = round(result3,3)
M = rbind(sl, result3)
rownames(M) = c("significance level","empirical type I error rate")
colnames(M) = 1:length(sl)
knitr::kable(M)

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
sl = c(0.01,0.05,0.1,0.2)
N = 3000 ##number of iterations
n = 30 ##sample size
result4 = rep(0,length(sl))

for(j in 1:length(sl)){
    
    t1e_vector = rep(0,n) 

    for(i in 1:N){
      x = runif(n, min = 0, max = 2)                                               
      q = qt(1-sl[j]/2, df = n-1)                                      
      t1e_vector[i] = mean(abs(mean(x)-1)>q*var(x)/sqrt(n))  ## the mean of chi^2(1) is 1
    }
    
    result4[j] = mean(t1e_vector)
  }

result4 = round(result4,3)
M = rbind(sl, result4)
rownames(M) = c("significance level","empirical type I error rate")
colnames(M) = 1:length(sl)
knitr::kable(M)

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
sl = c(0.01,0.05,0.1,0.2)
N = 3000 ##number of iterations
n = 30 ##sample size
result5 = rep(0,length(sl))

for(j in 1:length(sl)){
    
    t1e_vector = rep(0,n) 

    for(i in 1:N){
      x = rexp(n, rate = 1)                                               
      q = qt(1-sl[j]/2, df = n-1)                                      
      t1e_vector[i] = mean(abs(mean(x)-1)>q*var(x)/sqrt(n))  ## the mean of chi^2(1) is 1
    }
    
    result5[j] = mean(t1e_vector)
  }

result5 = round(result5,3)
M = rbind(sl, result5)
rownames(M) = c("significance level","empirical type I error rate")
colnames(M) = 1:length(sl)
knitr::kable(M)


## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
#PART1: plot
library(bootstrap)
par(mfrow = c(2, 5))

## ------------------------------------------------------------------------
#PART2: sample correlation matrix
print(cor(scor))

## ------------------------------------------------------------------------
#PART3: boostrap estimate
library(boot)
coreff <- function(x,i) cor(x[i,1],x[i,2])
M = matrix(0,4,1)

#standard errors of cor(mec,vec)
M[1,1] = sd(boot(data = scor[,c(1,2)], statistic = coreff, R=1000)$t)

#standard errors of cor(alg,ana)            
M[2,1] = sd(boot(data = scor[,c(3,4)], statistic = coreff, R=1000)$t)


#standard errors of cor(alg,sta)
M[3,1] = sd(boot(data = scor[,c(3,5)], statistic = coreff, R=1000)$t)


#standard errors of cor(ana,sta)
M[4,1] = sd(boot(data = scor[,c(4,5)], statistic = coreff, R=1000)$t)

colnames(M) = c("boot est of standard errors")
rownames(M) = c("cor(mec,vec)","cor(alg,ana)","cor(alg,sta)","cor(ana,sta)")
knitr::kable(M)

## ------------------------------------------------------------------------
b.sk <- function(x,i){
  m1 = mean(x[i,])
  m3 = mean((x[i,] - m1)^3) 
  m2 = mean((x[i,] - m1)^2) 
  return(m3 / m2^1.5) 
}


#PART1: Normal
N = 1000 ##number of iterations
n = 30   ##sample size
R = 1000 ##bootstrap number
cp.norm = rep(0,N)
cp.basic = rep(0,N)
cp.perc = rep(0,N)
count.norml = 0
count.normr = 0
count.basicl = 0
count.basicr = 0
count.percl = 0
count.percr = 0


for(i in 1:N){
  data_x = as.matrix(rnorm(n))
  obj = boot.ci(boot(data = data_x, statistic = b.sk, R = R), conf = 0.95, 
                type = c("norm", "basic", "perc"))
    
  cp.norm[i] = obj$norm[2] <= 0 & obj$norm[3] >= 0
  if (obj$norm[2] > 0) count.norml = count.norml+1
  if (obj$norm[3] < 0) count.normr = count.normr+1
  cp.basic[i] = obj$basic[4] <= 0 & obj$basic[5] >= 0
  if (obj$basic[4] > 0) count.basicl = count.basicl+1
  if (obj$basic[5] < 0) count.basicr = count.basicr+1
  cp.perc[i] = obj$perc[4] <= 0 & obj$perc[5] >= 0
  if (obj$perc[4] > 0) count.percl = count.percl+1
  if (obj$perc[5] < 0) count.percr = count.percr+1
  # the skewness of a normal distribution is 0 
  }
  
result1 = rbind(c(mean(cp.norm),mean(cp.basic),mean(cp.perc)),c(count.norml/N,count.basicl/N,count.percl/N),c(count.normr/N,count.basicr/N,count.percr/N))

#PART2: Chi^2(5)
N = 1000 ##number of iterations
n = 30   ##sample size
R = 1000 ##bootstrap number
cp.norm = rep(0,N)
cp.basic = rep(0,N)
cp.perc = rep(0,N)
count.norml = 0
count.normr = 0
count.basicl = 0
count.basicr = 0
count.percl = 0
count.percr = 0

for(i in 1:N){
  data_x = as.matrix(rchisq(n, df = 5))
  obj = boot.ci(boot(data = data_x, statistic = b.sk, R = R), conf = 0.95, 
                type = c("norm", "basic", "perc"))
    
  cp.norm[i] = obj$norm[2] <= 2*sqrt(0.4) & obj$norm[3] >= 2*sqrt(0.4)
  if (obj$norm[2] > 2*sqrt(0.4)) count.norml = count.norml+1
  if (obj$norm[3] < 2*sqrt(0.4)) count.normr = count.normr+1
  cp.basic[i] = obj$basic[4] <= 2*sqrt(0.4) & obj$basic[5] >= 2*sqrt(0.4)
  if (obj$basic[4] > 2*sqrt(0.4)) count.basicl = count.basicl+1
  if (obj$basic[5] < 2*sqrt(0.4)) count.basicr = count.basicr+1
  cp.perc[i] = obj$perc[4] <= 2*sqrt(0.4) & obj$perc[5] >= 2*sqrt(0.4)
  if (obj$perc[4] > 2*sqrt(0.4)) count.percl = count.percl+1
  if (obj$perc[5] < 2*sqrt(0.4)) count.percr = count.percr+1
  # the skewness of chi^2(fd) distribution is 2*sqrt(2/fd)
  }
  
result2 = rbind(c(mean(cp.norm),mean(cp.basic),mean(cp.perc)),c(count.norml/N,count.basicl/N,count.percl/N),c(count.normr/N,count.basicr/N,count.percr/N))

## ------------------------------------------------------------------------
M2 = rbind(result1,result2)
rownames(M2) = c("normal","proportion of miss left","proportion of miss right","chi^2(5)","proportion of miss left","proportion of miss right")
colnames(M2) = c("normal","basic","precentile")
knitr::kable(M2)

## ------------------------------------------------------------------------
pr1 = function (M){ ##the proportion of variance explained by the <U+FB01>rst principal component
  lambda = eigen(M)$values
  if (min(lambda)<=0) return(NULL)
  return(max(lambda)/sum(lambda))
}

jack_pr = function (data){
  n = ncol(data)
  theta.jack = rep(0,n)
  for (i in n){
    theta.jack[i] = pr1(M = cov(data[,-i]))
  }
  theta.hat = pr1(M = cov(data))
  return(c((n-1)*(mean(theta.jack)-theta.hat),sqrt((n-1)*mean((theta.jack - mean(theta.jack))^2))))
}

library(bootstrap)
result = jack_pr(scor)
M = matrix(data = result,1,2)
colnames(M) = c("est of bias","est of std")
rownames(M) = c("value")
knitr::kable(M)


## ------------------------------------------------------------------------
library(DAAG)
attach(ironslag) 

# for n-fold cross validation 
# fit models on leave-one-out samples 

jack_coef_cv = function(data = NULL){
  n = length(magnetic)
  e1 <- e2 <- e3 <- e4 <- rep(n,0)

for (i in 1:n) {
  y = magnetic[-i] 
  x = chemical[-i]
  
  J1 <- lm(y ~ x) 
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[i] 
  e1[i] <- magnetic[i] - yhat1
  
  J2 <- lm(y ~ x + I(x^2)) 
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[i] + J2$coef[3] * chemical[i]^2   
  e2[i] <- magnetic[i] - yhat2
  
  J3 <- lm(log(y) ~ x) 
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[i] 
  yhat3 <- exp(logyhat3)
  e3[i] <- magnetic[i] - yhat3
  
  J4 <- lm(y ~ poly(x,3,raw=T)) 
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[i] + J4$coef[3] * chemical[i]^2 + J4$coef[4] * chemical[i]^3
  e4[i] <- magnetic[i] - yhat4
}
  return(c(mean(e1^2),mean(e2^2),mean(e3^2),mean(e4^2)))
}

result = jack_coef_cv()
M = matrix(data = result,1,4)
colnames(M) = c("linear","quadratic","log","cubic")
rownames(M) = c("prediction error")
knitr::kable(M)

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
y = magnetic
x = chemical
summary(lm(y ~ x + I(x^2)) )

## ------------------------------------------------------------------------
R2_coef = function(data = NULL){
  e1 <- e2 <- e3 <- e4 <- 0

  y = magnetic
  x = chemical
  
  e1 <- summary(lm(y ~ x))$adj.r.squared
  e2 <- summary(lm(y ~ x + I(x^2)))$adj.r.squared
  e3 <- summary(lm(log(y) ~ x))$adj.r.squared
  e4 <- summary(lm(y ~ poly(x,3,raw=T)))$adj.r.squared

  return(c(e1,e2,e3,e4))
}

result = R2_coef()
M = matrix(data = result,1,4)
colnames(M) = c("linear","quadratic","log","cubic")
rownames(M) = c("Adjusted R Square")
knitr::kable(M)

## ------------------------------------------------------------------------

## ------------------------------------------------------------------------
y = magnetic
x = chemical
summary(lm(log(y) ~ x))

## ------------------------------------------------------------------------
maxout = function(x,y){
  X = x - mean(x) 
  Y = y - mean(y) 
  outx = sum(X > max(Y)) + sum(X < min(Y)) 
  outy = sum(Y > max(X)) + sum(Y < min(X)) 
  return(max(c(outx, outy)))
}
count5_permutation = function(X,Y,N){
  if(length(X)<=length(Y)){
    m = length(X)
    Z = c(X,Y)
    theta.hat = maxout(X,Y)
  }
  else {
    m = length(Y)
    Z = c(Y,X)
    theta.hat = maxout(Y,X)
  }
  theta.per = rep(0,N)
  for (i in 1:N){
    ind = sample(1:length(Z),size = 2*m, replace = F)
    z0 = Z[ind]
    theta.per[i] = maxout(z0[1:m],z0[(m+1):(2*m)])
  }
  out = (1+sum(theta.per>=theta.hat))/(1+N)  ##p-value of count-five test
}

## ------------------------------------------------------------------------
##exmaple 6.15

m = 500
n = 600
N0 = 300
result = rep(0,N0)
for(i in 1:N0){
  result[i] = count5_permutation(rnorm(m),rnorm(n),1e3)
}


## ------------------------------------------------------------------------
mean(result) ##p-value

## ------------------------------------------------------------------------
##exmaple 6.15

m = 500
n = 600
N0 = 300
result = rep(0,N0)
for(i in 1:N0){
  result[i] = count5_permutation(rnorm(m),1.5*rnorm(n),1e3)
}


## ------------------------------------------------------------------------
mean(result) ##p-value

## ------------------------------------------------------------------------
library(Ball)
library(energy)

power_comparison = function(n, N, alpha = 0.1){
  ## n is the sample size
  ## N is the number of replication
  p.dist = matrix(0,N,2)
  p.ball = matrix(0,N,2)
  
  for(i in 1:N){
    X = matrix(rnorm(2*n), nrow = n, ncol = 2)
    e = matrix(rnorm(2*n), nrow = n, ncol = 2)
    Y1 = X/4 + e
    Y2 = X/4 * e
    p.ball[i,1] = bd.test(x = X, y = Y1, seed = i*54321, R = 99)$p.value
    p.dist[i,1] = dcor.test(x = X, y = Y1, R = 99)$p.value
    p.ball[i,2] = bd.test(x = X, y = Y2, seed = i*12345, R = 999)$p.value
    p.dist[i,2] = dcor.test(x = X, y = Y2, R = 99)$p.value
  }
  return(c(mean(p.ball[,1] < alpha), mean(p.dist[,1] < alpha), mean(p.ball[,2] < alpha), mean(p.dist[,2] < alpha)))
}


vec.n = seq(10, 60, 10)
matrix.power = matrix(0, nrow = 5, ncol = length(vec.n))
for(i in 1:length(vec.n)){
  matrix.power[,i] = c(vec.n[i],power_comparison(vec.n[i], 2e2))
}
# rownames(matrix.power) = c( "sample size","dcor test1", "ball test1","dcor test2", "ball test2" )
# knitr::kable(matrix.power)

plot(vec.n, matrix.power[2,], type="l", col = "red", lty=3, ylim = c(0,1),xlab = "sample size", ylab = "power", main = "model1")
lines(vec.n, matrix.power[3,], type = "l", col = "blue", lty =3)
legend("bottomright", c("ball","dcor"), fill = c("red","blue"))

plot(vec.n, matrix.power[4,], type="l", col = "red", lty=3, ylim = c(0,1),xlab = "sample size", ylab = "power", main = "model2")
lines(vec.n, matrix.power[5,], type = "l", col = "blue", lty =3)
legend("bottomright", c("ball","dcor"), fill = c("red","blue"))

## ------------------------------------------------------------------------
v = c(0.1,1,10,100)
M = matrix(0,2,length(v))
for (i in 1:length(v)){
  M[1,i] = identical(log(exp(v[i])),exp(log(v[i])))
  M[2,i] = isTRUE(all.equal(log(exp(v[i])),exp(log(v[i]))))
}

## ----echo=FALSE----------------------------------------------------------
colnames(M) = v
rownames(M) = c("identical","all.equal")
knitr::kable(M)

## ------------------------------------------------------------------------
func_S = function(a,k){
  b = 1-pt(q = sqrt((a^2)*k/(k+1-a^2)),df = k)
  return(b)
}
vec_k = c(5:10,50:55)
result2 = rep(0,length(vec_k))
for(i in 1:length(vec_k)){
  result2[i] = uniroot(function(y){func_S(a=y,k=vec_k[i])-func_S(a=y,k=vec_k[i]-1)},interval = c(1,min(2,sqrt(vec_k[i]))))$root
}

## ------------------------------------------------------------------------
func_intel = function(a,k){
  f = function(u) (1+u^2/k)^(-(k+1)/2)
  b = exp(lgamma((k+1)/2)-lgamma(k/2))/sqrt(k)*integrate(f,lower = 0,upper = sqrt((a^2)*k/(k+1-a^2)))$value
  return(b)
}

vec_k = c(5:10,50:55)
result = rep(0,length(vec_k))
for(i in 1:length(vec_k) ){
  result[i] = uniroot(function(y) {func_intel(a = y,k = vec_k[i])-func_intel(a = y,k = vec_k[i]-1)},interval = c(1,min(2,sqrt(vec_k[i])-1e-4)))$root
}

## ----echo=FALSE----------------------------------------------------------
M2 = matrix(0,2,length(vec_k))
M2[1,] = result
M2[2,] = result2
colnames(M2) = vec_k
rownames(M2) = c("Ex5","Ex4")
knitr::kable(M2)

## ------------------------------------------------------------------------
E_step = function(para,p0,q0,nA=28,nB=24,nO=41,nAB=70){
  p = para[1]
  q = para[2]
  b = 2*nA*log(p)+2*nB*log(q)+2*nO*log(1-p-q)+nAB*log(p)+nAB*log(q)+nA*2*(1-p0-q0)/(2-2*q0-p0)*log(2*(1-p-q)/p)+nB*2*(1-p0-q0)/(2-2*p0-q0)*log(2*(1-p-q)/q)
  return(b)
}

##M_step
N = 10
p0 = 0.2
q0 = 0.2
p_rec = rep(0,N)
q_rec = rep(0,N)
loglike = rep(0,N)

for(i in 1:N){
  obj = function(y) -E_step(para = y, p0,q0)
  itr = optim(par = c(0.3,0.3),fn = obj)
  p0 = itr$par[1]
  q0 = itr$par[2]
  p_rec[i] = itr$par[1]
  q_rec[i] = itr$par[2]
  loglike[i] = itr$value
}

## ----echo=FALSE----------------------------------------------------------
plot(1:N,loglike,col = "red",type = "l",xlab = "iteration",ylab = "log-likelihood")

plot(p_rec,col = "blue",type = "l",ylim = c(0.28,0.35),xlab = "iteration",ylab = "p&q")
lines(q_rec,col = "purple",type = "l")
legend("topright",c("p","q"),lty = c(1,1),col = c("blue","purple"))
M3 = matrix(0,1,2)
M3[1,] = round(c(p_rec[N],q_rec[N]),3)
colnames(M3) = c("p","q")
rownames(M3) = c("MLE")
knitr::kable(M3)


## ------------------------------------------------------------------------
formulas <- list( 
  mpg ~ disp, 
  mpg ~ I(1 / disp), 
  mpg ~ disp + wt, 
  mpg ~ I(1 / disp) + wt ) 

## ------------------------------------------------------------------------
## for loop
for(i in 1:length(formulas)){
  result = lm(formulas[[i]], data = mtcars)
  print(result)
}

## ------------------------------------------------------------------------
## lapply
lapply(formulas, function(f) lm(f, data = mtcars))

