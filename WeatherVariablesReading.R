#setwd("C:/Users/amalonso/AMAlonso-1/Alumnos/AA Zehang Li/TFM_method_in_day-ahead_market/ExogenousVariables")
setwd("D:/Script/R_script/supply_curves_functional_prediction/ExogenousVariables")
Files = dir()

library(readr)
WeatherVariables = read_csv(Files[1], skip = 4, col_types = cols(X1 = col_skip()), col_names = FALSE)

for (i in 2:52){
  WeatherVariables = cbind(WeatherVariables, 
                           read_csv(Files[i], col_names = FALSE, 
                                    col_types = cols(X1 = col_skip()), skip = 4))
}

library(stats)
library(MASS)

# Explain variability.
for (i in 1:8){
  PCA = summary(princomp(WeatherVariables[,seq(i,416,8)]))
  EV = cumsum(PCA$sdev^2/sum(PCA$sdev^2))
  if (i == 1){
    n.pc = which(EV >= 0.9)[1]
  } 
  else {
    n.pc = c(n.pc, which(EV >= 0.9)[1])
  }
}
n.pc

# Selected number of principal components. 
n.pc = c(3, 0, 9, 9, 9, 3, 3, 0)

for (i in c(1, 3:7)){
  PCA = princomp(WeatherVariables[,seq(i,416,8)])
  if (i == 1){
    Weather.pc = as.matrix(PCA$scores[,1:n.pc[i]])
  } 
  else {
    Weather.pc = cbind(Weather.pc, as.matrix(PCA$scores[,1:n.pc[i]]))
  }
}
dim(Weather.pc)

save(Weather.pc, file = "Weather_pc.RData")

for (i in c(1, 3:6)){
  WM = apply(WeatherVariables[,seq(i,416,8)], 1, mean)
  WS = apply(WeatherVariables[,seq(i,416,8)], 1, sd)
  if (i == 1){
    Weather.stat = cbind(WM, WS)
  } 
  else {
    Weather.stat = cbind(Weather.stat, WM, WS)
  }
}
dim(Weather.stat)

save(Weather.stat, file = "Weather_Stat.RData")


