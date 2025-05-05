library(vroom)
library(fda)
library(fda.usc)
library(foreach)
library(forecast)
library(doParallel)
library(doSNOW)
library(tidyverse)
library(ggplot2)
library(tictoc)

#################################################################################
# 1. loading data, evaluating curves at grids and representing by order-1 B-splines 
#################################################################################

# load("D:/Script/R_script/day_ahead market/day_ahead market/day_ahead_supply_price.RData")
# load("D:/Script/R_script/day_ahead market/day_ahead market/day_ahead_supply_cumsum.RData")

load("C:/Users/amalonso/AMAlonso-1/Alumnos/AA Zehang Li/DayAheadMarketDistanceLearning/day_ahead_supply_price.RData")
load("C:/Users/amalonso/AMAlonso-1/Alumnos/AA Zehang Li/DayAheadMarketDistanceLearning/day_ahead_supply_cumsum.RData")

Prices = day_ahead_supply_price[1:(365*24*5 +48),3:700]
Quantities = day_ahead_supply_cumsum[1:(365*24*5 +48),3:700]
rm(list=c('day_ahead_supply_price','day_ahead_supply_cumsum'))

myapprox=function(x,y,evals){
  fc <- approxfun(x, y, method = "const",f = 1,rule=2)
  result=fc(evals)
}

regular.nodes = seq(0, 180.3, length = 200)
Qapprox=sapply(1:nrow(Quantities), function(i) myapprox(Prices[i,],Quantities[i,],regular.nodes))

rng=c(min(regular.nodes),max(regular.nodes))
basis=create.bspline.basis(rng,nbasis = 20,norder = 2)
quantities.fd=Data2fd(regular.nodes,Qapprox,basis)
plot(quantities.fd[30], 
     ylim = c(min(Qapprox[,30])-1000,max(Qapprox[,30])+1000))
lines(regular.nodes,Qapprox[,30],col="red",type = "s")


Testset = Qapprox[,-(1:(365 *24 *4 +24))]


##############################################################
# 2. functional regression 
##############################################################


acfs=sapply(1:nrow(Qapprox), function(i) acf(Qapprox[i,],lag.max = 24*7,plot = F)$acf)

acfs.df=as.data.frame(acfs)
acfs.df=gather(acfs.df,everything(),key="lines",value="value")
acfs.df$price=rep(regular.nodes,each=169)
acfs.df$lags=rep(0:168,200)


gg0=ggplot(data=acfs.df,mapping=aes(x=lags,y=value))+
  geom_line(aes(group=lines,color=price))+scale_color_distiller(trans="reverse")+
  labs(y="ACF",color = "Quantity")+
  scale_x_continuous(breaks=seq(0, 170, 5))+
  scale_y_continuous(breaks=seq(-1, 1, 0.25))+
  labs(title = toupper("ACF's of time series fixed at 500 quantities"))+
  theme_bw() + theme(legend.position = "right",plot.title = element_text(hjust = 0.5),
                     panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_vline(xintercept = 3) +
  geom_vline(xintercept = 22) +
  geom_vline(xintercept = 26);gg0


tic()

start = 0
end = 365

cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = (end - start +1) * 24
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

quantity_FunReg_Noassump= foreach(i=start:end,.combine = 'cbind',.packages = c("fda","fda.usc"),.options.snow = opts,.errorhandling="pass")%dopar%{
  #foreach(j=1:24,.combine = "cbind",.packages = c("fda","fda.usc"),.options.snow = opts,.errorhandling="pass")%dopar%{
  yhat = matrix(data = NA, nrow = 200, ncol = 24)

  Time = 24*365*4 + 24 + 24*i 
  Win_size = 24*30*4

  # The fRegress fails when norder = 1 even in the simplest case, y ~ x1.
  XXbasis=create.bspline.basis(rng,nbasis = 20,norder = 2)
  XX1 <- Qapprox[,(Time-Win_size+1):Time]
  Lags = c(1,2,3,23,24,25)
  indexes = seq(max(Lags)+1,Win_size,1) # for y
  y = XX1[,indexes]
  y= smooth.basisPar(regular.nodes,y, XXbasis)$fd
  k = 1
  for (l in Lags){
    xj = XX1[,indexes-l]
    xj = smooth.basisPar(regular.nodes,xj, XXbasis)$fd
    eval(str2expression(paste0("x", k, "= xj")))
    eval(str2expression(paste0("beta", k, "= with(x", k, ", fd(basisobj=basis, fdnames=fdnames))")))
    k = k + 1
  }
  xfdlist  = list(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6)
  betalist = list(x1=fdPar(beta1), x2=fdPar(beta2), x3=fdPar(beta3), 
                  x4=fdPar(beta4), x5=fdPar(beta5), x6=fdPar(beta6))
      
  fRegressout <- fRegress(y, xfdlist, betalist, method = "fRegress")
  for (j in 1:24){
    k = 1
    for (l in Lags){
      xj = XX1[,max(indexes)-l+j]
      xj = smooth.basisPar(regular.nodes,xj, XXbasis)$fd
      eval(str2expression(paste0("x", k, "= xj")))
      k = k + 1
    }
    xfdlist  = list(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6)
    yhat0 = predict.fRegress(fRegressout, xfdlist)
    yhat[,j]= eval.fd(yhat0,regular.nodes)
    XX1 = cbind(XX1, yhat[,j])
  }
  yhat = yhat
}
toc()
close(pb)
stopCluster(cl) 

save(quantity_FunReg_Noassump, file = 'quantity_FunReg_recursive.RData')

plot(regular.nodes,quantity_FunReg_Noassump[,96], type = "s", 
     ylim = c(min(Qapprox[,(24*365*4 + 24  +96)])-1000, 
              max(Qapprox[,(24*365*4 + 24  +96)])+1000))
lines(regular.nodes,Qapprox[,(24*365*4 + 24  +96)],col="red",type="s")


quantity_FunReg_Noassump = do.call(cbind,mget(ls()[startsWith(ls(),"quantity_FunReg_Noassump")]))
save(quantity_FunReg_Noassump, file = "D:/Script/R_script/supply_curves_functional_prediction/quantity_FunReg_Noassump.RData")
#rm(list = ls()[startsWith(ls(),"quantity_FunReg_Noassump")][2:14])

##############################################################
# 3. PCA + ARIMA 
##############################################################

tic()
start = 0
end = 365 

cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = (end - start +1) * 24
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

quantity_PrePCArima = foreach(i=start:end,.combine = "cbind",.packages = c("fda","fda.usc","forecast"),.options.snow = opts)%dopar%{
  yhat = matrix(data = NA, nrow = 200, ncol = 24)
  #foreach(j=1:24,.combine = "cbind",.packages = c("fda","fda.usc","forecast"),.options.snow = opts)%dopar%{
  
  Time = 24*365*4 + 24 + 24*i 
  width = 24*30
  
  fpca <- pca.fd(fdobj=quantities.fd[(Time-width-1):Time], nharm = 10, centerfns=TRUE)
  scores=fpca$scores
  weight_functions <- fpca$harmonics
  mean_functions = fpca$meanfd
  
  ts1 <- ts(scores[,1],frequency = 24)  
  ts2 <- ts(scores[,2],frequency = 24)
  ts3 <- ts(scores[,3],frequency = 24)
  ts4 <- ts(scores[,4],frequency = 24)
  ts5 <- ts(scores[,5],frequency = 24)
  ts6 <- ts(scores[,6],frequency = 24)
  ts7 <- ts(scores[,7],frequency = 24)
  ts8 <- ts(scores[,8],frequency = 24)
  ts9 <- ts(scores[,9],frequency = 24)
  ts10 <- ts(scores[,10],frequency = 24)
  
  fit1=auto.arima(ts1,d=0)
  fit2=auto.arima(ts2,d=0)
  fit3=auto.arima(ts3,d=0)
  fit4=auto.arima(ts4,d=0)
  fit5=auto.arima(ts5,d=0)
  fit6=auto.arima(ts6,d=0)
  fit7=auto.arima(ts7,d=0)
  fit8=auto.arima(ts8,d=0)
  fit9=auto.arima(ts9,d=0)
  fit10=auto.arima(ts10,d=0)

  for (j in 1:24){  
    K = j #K step ahead
    scores.pre1=as.numeric(forecast(fit1,h=K)$mean)
    scores.pre2=as.numeric(forecast(fit2,h=K)$mean)
    scores.pre3=as.numeric(forecast(fit3,h=K)$mean)
    scores.pre4=as.numeric(forecast(fit4,h=K)$mean)
    scores.pre5=as.numeric(forecast(fit5,h=K)$mean)
    scores.pre6=as.numeric(forecast(fit6,h=K)$mean)
    scores.pre7=as.numeric(forecast(fit7,h=K)$mean)
    scores.pre8=as.numeric(forecast(fit8,h=K)$mean)
    scores.pre9=as.numeric(forecast(fit9,h=K)$mean)
    scores.pre10=as.numeric(forecast(fit10,h=K)$mean)

    a1 <- scores.pre1[K]*weight_functions[1]
    a2 <- scores.pre2[K]*weight_functions[2]
    a3 <- scores.pre3[K]*weight_functions[3]
    a4 <- scores.pre4[K]*weight_functions[4]
    a5 <- scores.pre5[K]*weight_functions[5]
    a6 <- scores.pre6[K]*weight_functions[6]
    a7 <- scores.pre7[K]*weight_functions[7]
    a8 <- scores.pre8[K]*weight_functions[8]
    a9 <- scores.pre9[K]*weight_functions[9]
    a10 <- scores.pre10[K]*weight_functions[10]
    
    at <- a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + mean_functions
    yhat[,j] = eval.fd(regular.nodes,at)
  }
  return(yhat)
}
close(pb)
stopCluster(cl) 
toc()


plot(regular.nodes,quantity_PrePCArima[,3])
lines(regular.nodes,Qapprox[,(24*365*4 + 24  +3)],col="red")


#############################################################################
# 4. stepwisezation, replacing the prediction with the closest curve in history
#############################################################################


DistanceCalculation = function(Price1,Quant1,Price2,Quant2,p){
  Price1 = c(0, Price1) # The point (0,0) is necessary for the distance calculation.
  Quant1 = c(0, Quant1)

  Price2 = c(0, Price2)
  Quant2 = c(0, Quant2)
  
  # This part implements the elimination of redundant offers with 
  # the same price. We take the maximum cumulative quantity of 
  # the offers with the same price. With the new data, this step 
  # isn't necessary.
  UPrice1 = intersect(Price1,Price1)
  # UQuant1 = UPrice1*NA
  # for (i in 1:length(UPrice1)){
  #   UQuant1[i] = max(Quant1[Price1==UPrice1[i]])
  # }
  UPrice2 = intersect(Price2,Price2)
  # UQuant2 = UPrice2*NA
  # for (i in 1:length(UPrice2)){
  #   UQuant2[i] = max(Quant2[Price2==UPrice2[i]])
  # }
  UPrice1 = Price1[-1]
  UPrice2 = Price2[-1]
  UQuant1 = Quant1[-1]
  UQuant2 = Quant2[-1]
  
  # The next lines are the distance calculations. The previous should 
  # be calculated in advance. The input of the distance calculation 
  # are the vectors UQuant1, UQuant2, UPrice1 and UPrice2.
  
  # Where the jumps of the two functions appears:
  T = sort(c(UQuant1,UQuant2)) # This can be optimized.
  N = length(T)
  T = T[2:(N-1)]
  N = N-2
  
  # Auxiliary functions:
  Orders <- function(c, x){
    O <- sum(c >= x)
    return(O)
  } 
  # The area of a rectangle equals base x height.
  Bases = diff(T)
  Centers = as.matrix((T[2:N] + T[1:(N-1)])/2)
  O1 = apply(Centers, 1, Orders, UQuant1) 
  O2 = apply(Centers, 1, Orders, UQuant2) 
  #p = 1 # p = 1 para L1, p = 2 para L2, ...
  
  H1 = UPrice1[O1 + (O1 < length(UPrice1))] # when O1+1 is larger than the length, eliminating the +1
  H2 = UPrice2[O2 + (O2 < length(UPrice2))]
  
  Heights = abs(H1 - H2)^p
  Distance = sum(Bases*Heights)
  return(Distance)
}

PQ.DistanceCalculation = function(Price1,Quant1,Price2,Quant2,p){
  # This distance assume that Price1 = Price2 = regular.nodes
  Distance = sum(abs(Quant1 - Quant2)^p)*Price1[2] # Price1[2] = dp
  if (p == 2){Distance = sqrt(Distance)}
  return(Distance)
}  

closest.curve=function(curve,t){
  res1 = t %% 24
  res2 = which(1:(365 *24 *4 +24 + (t-24)) %% 24 == res1)
  DM = sapply(res2,function(i) PQ.DistanceCalculation(regular.nodes,curve,regular.nodes,Qapprox[,i],p=1)) 
  index=res2[which.min(DM)]
  result=Qapprox[,index]
  return(result)
}

PQ.closest.curve=function(curve,t){
  res1 = t %% 24
  res2 = 1:(365 *24 *4 + 24 + t-24)#which(1:(365 *24 *4 +24) %% 24 == res1)
  DM = sapply(res2,function(i) PQ.DistanceCalculation(regular.nodes,curve,regular.nodes,Qapprox[,i],p=1)) 
  index=res2[which.min(DM)]
  result=Qapprox[,index]
  return(result)
}

# Finding the closest curve for functional regression.
tic()
cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = ncol(quantity_FunReg_Noassump)
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

quantity_FunReg_Noassump_closest  = foreach(i = 1:iterations, .combine = "cbind", .options.snow = opts)%dopar%{
  PQ.closest.curve(quantity_FunReg_Noassump[,i],i)}
close(pb)
stopCluster(cl) 
toc()

plot(regular.nodes,quantity_FunReg_Noassump[,3])
lines(regular.nodes,Qapprox[,(24*365*4 + 24  +3)],col="red")
lines(regular.nodes,quantity_FunReg_Noassump_closest[,3],col='blue')

# Finding the closest curve for PCA+ARIMA.
tic()
cl = makeSOCKcluster(6)
registerDoSNOW(cl)
iterations = ncol(quantity_FunReg_Noassump)
pb = txtProgressBar(max=iterations, style=3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress=progress)

quantity_PrePCArima_Noassump_closest  = foreach(i = 1:iterations, .combine = "cbind", .options.snow = opts)%dopar%{
  PQ.closest.curve(quantity_PrePCArima_Noassump[,i],i)}
close(pb)
stopCluster(cl) 
toc()

plot(regular.nodes,quantity_PrePCArima_Noassump[,3])
lines(regular.nodes,Qapprox[,(24*365*4 + 24  +3)],col="red")
lines(regular.nodes,quantity_PrePCArima_Noassump_closest[,3],col='blue')



#############################################################################
# 5.Comparing with Naive prediction
#############################################################################

result.naive=cbind(Qapprox[,1:24],Qapprox[,1:(ncol(Qapprox)-24)])
result.naive = result.naive[,-(1:(365 *24 *4 +24))]

err_naive = sapply(1:ncol(result.naive),function(i) 
  {PQ.DistanceCalculation(regular.nodes,result.naive[,i],regular.nodes,Testset[,i],p=1)})
mean(err_naive)
sd(err_naive)
median(err_naive)
mad(err_naive)

err_FunReg = sapply(1:ncol(quantity_FunReg_Noassump),function(i) 
  {PQ.DistanceCalculation(regular.nodes,quantity_FunReg_Noassump[,i],regular.nodes,Testset[,i],p=1)})
mean(err_FunReg)
sd(err_FunReg)
median(err_FunReg)
mad(err_FunReg)

err_closest1 = sapply(1:ncol(quantity_FunReg_Noassump_closest),function(i) 
{PQ.DistanceCalculation(regular.nodes,quantity_FunReg_Noassump_closest[,i],regular.nodes,Testset[,i],p=1)})
mean(err_closest1)
sd(err_closest1)
median(err_closest1)
mad(err_closest1)

order(err_FunReg,decreasing = T)
plot(regular.nodes,Testset[,123],ylim =c(25000,100000))
lines(regular.nodes,quantity_FunReg_Noassump[,123],col="red")
lines(regular.nodes,result.naive[,123],col='blue')
#5883 - 5902 are not monotonic


load("D:/Script/R_script/supply_curves_functional_prediction/quantity_PrePCArima_Noassump.RData")
err_PCArima = sapply(1:ncol(quantity_PrePCArima_Noassump),function(i) {PQ.DistanceCalculation(regular.nodes,quantity_PrePCArima_Noassump[,i],regular.nodes,Testset[,i],p=1)})
mean(err_PCArima)
sd(err_PCArima)
median(err_PCArima)
mad(err_PCArima)

order(err_PCArima,decreasing = T)
plot(regular.nodes,Testset[,123],ylim =c(25000,100000))
lines(regular.nodes,quantity_PrePCArima_Noassump[,123],col="red")
lines(regular.nodes,result.naive[,123],col='blue')

load("D:/Script/R_script/supply_curves_functional_prediction/quantity_PrePCArima_Noassump_closest.RData")

err_closest2 = sapply(1:ncol(quantity_PrePCArima_Noassump_closest),function(i) {PQ.DistanceCalculation(regular.nodes,quantity_PrePCArima_Noassump_closest[,i],regular.nodes,Testset[,i],p=1)})
mean(err_closest2)
sd(err_closest2)
median(err_closest2)
mad(err_closest2)
