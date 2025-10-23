rm(list=ls())

library(iterators)
library(parallel)
library(foreach)
library(doParallel)
library(MCS)

cl <- makeCluster(14)
registerDoParallel(cl)

setwd("E:/software/R/code/CoVaRCoes")
result <- matrix(nrow=31,ncol=9)
MCSresult95M<- matrix(nrow=1,ncol=9)
MCSresult90M<- matrix(nrow=1,ncol=9)
MCSresult95R<- matrix(nrow=1,ncol=9)
MCSresult90R<- matrix(nrow=1,ncol=9)



test <- read.csv("E:/software/R/code/CoVaRCoes/Gaui2sfullCoVaR.csv",head=TRUE)
data <- test 
system.time(result <- foreach(i=1:31,.packages = "MCS")%dopar%{
  l <- 1
  set.seed(100*i)
  datai <- data[,(9*(i-1)+1):(9*i)]
  MCS <- MCSprocedure(Loss = datai, alpha = 0, B = 5000, statistic = "Tmax")
  a<- MCS@show[,3]
  b<- MCS@show[,6]
  while (length(b)<9){
    set.seed(100*i+100*l)
    MCS <- MCSprocedure(Loss = datai, alpha = -0.001, B = 10000, statistic = "Tmax")
    a<- MCS@show[,3]
    b<- MCS@show[,6]
    l <- l+1
  }
  print(cbind(a,b))
})

MCSM <- matrix(nrow=9,ncol=31)
MCSR <- matrix(nrow=9,ncol=31)
for (j in 1:31){
  r<- result[[j]]
  MCSM[,j] <- r[,1]
  MCSR[,j] <- r[,2]
}

for (i in 1:9){
  MCSresult90M[i] <- length(MCSM[i,MCSM[i,]>0.90])
  MCSresult95M[i] <- length(MCSM[i,MCSM[i,]>0.95])
  MCSresult90R[i] <- length(MCSR[i,MCSR[i,]>0.90])
  MCSresult95R[i] <- length(MCSR[i,MCSR[i,]>0.95])
}
Ai2s_89=t(rbind(MCSresult95R,MCSresult90R))




test <- read.csv("E:/software/R/code/CoVaRCoes/Gaus2ifullCoVaR.csv",head=TRUE)
data <- test 
system.time(result <- foreach(i=1:31,.packages = "MCS")%dopar%{
  l <- 1
  set.seed(100*i)
  datai <- data[,(9*(i-1)+1):(9*i)]
  MCS <- MCSprocedure(Loss = datai, alpha = 0, B = 5000, statistic = "Tmax")
  a<- MCS@show[,3]
  b<- MCS@show[,6]
  while (length(b)<9){
    set.seed(100*i+100*l)
    MCS <- MCSprocedure(Loss = datai, alpha = -0.001, B = 10000, statistic = "Tmax")
    a<- MCS@show[,3]
    b<- MCS@show[,6]
    l <- l+1
  }
  print(cbind(a,b))
})

MCSM <- matrix(nrow=9,ncol=31)
MCSR <- matrix(nrow=9,ncol=31)
for (j in 1:31){
  r<- result[[j]]
  MCSM[,j] <- r[,1]
  MCSR[,j] <- r[,2]
}

for (i in 1:9){
  MCSresult90M[i] <- length(MCSM[i,MCSM[i,]>0.90])
  MCSresult95M[i] <- length(MCSM[i,MCSM[i,]>0.95])
  MCSresult90R[i] <- length(MCSR[i,MCSR[i,]>0.90])
  MCSresult95R[i] <- length(MCSR[i,MCSR[i,]>0.95])
}
As2i_89=t(rbind(MCSresult95R,MCSresult90R))



test <- read.csv("E:/software/R/code/CoVaRCoes/POTi2sfullCoVaR.csv",head=TRUE)
data <- test 
system.time(result <- foreach(i=1:31,.packages = "MCS")%dopar%{
  l <- 1
  set.seed(100*i)
  datai <- data[,(9*(i-1)+1):(9*i)]
  MCS <- MCSprocedure(Loss = datai, alpha = 0, B = 5000, statistic = "Tmax")
  a<- MCS@show[,3]
  b<- MCS@show[,6]
  while (length(b)<9){
    set.seed(100*i+100*l)
    MCS <- MCSprocedure(Loss = datai, alpha = -0.001, B = 10000, statistic = "Tmax")
    a<- MCS@show[,3]
    b<- MCS@show[,6]
    l <- l+1
  }
  print(cbind(a,b))
})

MCSM <- matrix(nrow=9,ncol=31)
MCSR <- matrix(nrow=9,ncol=31)
for (j in 1:31){
  r<- result[[j]]
  MCSM[,j] <- r[,1]
  MCSR[,j] <- r[,2]
}

for (i in 1:9){
  MCSresult90M[i] <- length(MCSM[i,MCSM[i,]>0.90])
  MCSresult95M[i] <- length(MCSM[i,MCSM[i,]>0.95])
  MCSresult90R[i] <- length(MCSR[i,MCSR[i,]>0.90])
  MCSresult95R[i] <- length(MCSR[i,MCSR[i,]>0.95])
}
Bi2s_89=t(rbind(MCSresult95R,MCSresult90R))


test <- read.csv("E:/software/R/code/CoVaRCoes/POTs2ifullCoVaR.csv",head=TRUE)
data <- test 
system.time(result <- foreach(i=1:31,.packages = "MCS")%dopar%{
  l <- 1
  set.seed(100*i)
  datai <- data[,(9*(i-1)+1):(9*i)]
  MCS <- MCSprocedure(Loss = datai, alpha = 0, B = 5000, statistic = "Tmax")
  a<- MCS@show[,3]
  b<- MCS@show[,6]
  while (length(b)<9){
    set.seed(100*i+100*l)
    MCS <- MCSprocedure(Loss = datai, alpha = -0.001, B = 10000, statistic = "Tmax")
    a<- MCS@show[,3]
    b<- MCS@show[,6]
    l <- l+1
  }
  print(cbind(a,b))
})

MCSM <- matrix(nrow=9,ncol=31)
MCSR <- matrix(nrow=9,ncol=31)
for (j in 1:31){
  r<- result[[j]]
  MCSM[,j] <- r[,1]
  MCSR[,j] <- r[,2]
}

for (i in 1:9){
  MCSresult90M[i] <- length(MCSM[i,MCSM[i,]>0.90])
  MCSresult95M[i] <- length(MCSM[i,MCSM[i,]>0.95])
  MCSresult90R[i] <- length(MCSR[i,MCSR[i,]>0.90])
  MCSresult95R[i] <- length(MCSR[i,MCSR[i,]>0.95])
}
Bs2i_89=t(rbind(MCSresult95R,MCSresult90R))

##----------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------
##
##                            Table 3 colume 7-8
##          Panel A Industries → System colume 7-8
Ai2s_89

##          Panel A System → Industries colume 7-8
As2i_89

##          Panel B Industries → System colume 7-8
Bi2s_89

##          Panel B System → Industries colume 7-8
Bs2i_89



