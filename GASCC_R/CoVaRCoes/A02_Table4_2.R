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

MCSresult95R<- matrix(nrow=1,ncol=9)
MCSresult90R<- matrix(nrow=1,ncol=9)

## --------------------------------------------------------------------------------------------------------------------------------
##          Panel A  Gaussian   Industries → System
PAi2s <- matrix(nrow=9,ncol=9)
testall <- read.csv("E:/software/R/code/CoVaRCoes/Gaui2sfullCoVaRCoES.csv",head=TRUE)
for (k in 1:3){
  data <- testall[,(279*(k-1)+1):(279*k)]
result <- foreach(i=1:31,.packages = "MCS")%dopar%{
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
}

MCSR <- matrix(nrow=9,ncol=31)
for (j in 1:31){
  r<- result[[j]]
  MCSR[,j] <- r[,2]
}

for (i in 1:9){
  MCSresult90R[i] <- length(MCSR[i,MCSR[i,]>0.90])
  MCSresult95R[i] <- length(MCSR[i,MCSR[i,]>0.95])
}
A=t(rbind(MCSresult95R,MCSresult90R))


PAi2s[,(3*(k-1)+1):(3*k)] <- cbind(t(t(10*rowMeans(matrix(as.matrix(colMeans(as.matrix(data))), nrow = 9, ncol = 31)))),A)
}

## --------------------------------------------------------------------------------------------------------------------------------
## Panel A  Gaussian    System → Industries
PAs2i <- matrix(nrow=9,ncol=9)
testall <- read.csv("E:/software/R/code/CoVaRCoes/Gaus2ifullCoVaRCoES.csv",head=TRUE)
for (k in 1:3){
  data <- testall[,(279*(k-1)+1):(279*k)]
  result <- foreach(i=1:31,.packages = "MCS")%dopar%{
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
  }
  
  MCSR <- matrix(nrow=9,ncol=31)
  for (j in 1:31){
    r<- result[[j]]
    MCSR[,j] <- r[,2]
  }
  
  for (i in 1:9){
    MCSresult90R[i] <- length(MCSR[i,MCSR[i,]>0.90])
    MCSresult95R[i] <- length(MCSR[i,MCSR[i,]>0.95])
  }
  A=t(rbind(MCSresult95R,MCSresult90R))
  
  PAs2i[,(3*(k-1)+1):(3*k)] <- cbind(t(t(10*rowMeans(matrix(as.matrix(colMeans(as.matrix(data))), nrow = 9, ncol = 31)))),A)
}

##===================================================================================================
##          Panel B  POT   Industries → System
PBi2s <- matrix(nrow=9,ncol=9)
testall <- read.csv("E:/software/R/code/CoVaRCoes/POTi2sfullCoVaRCoES.csv",head=TRUE)
for (k in 1:3){
  data <- testall[,(279*(k-1)+1):(279*k)]
  result <- foreach(i=1:31,.packages = "MCS")%dopar%{
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
  }
  
  MCSR <- matrix(nrow=9,ncol=31)
  for (j in 1:31){
    r<- result[[j]]
    MCSR[,j] <- r[,2]
  }
  
  for (i in 1:9){
    MCSresult90R[i] <- length(MCSR[i,MCSR[i,]>0.90])
    MCSresult95R[i] <- length(MCSR[i,MCSR[i,]>0.95])
  }
  A=t(rbind(MCSresult95R,MCSresult90R))
  
  PBi2s[,(3*(k-1)+1):(3*k)] <- cbind(t(t(10*rowMeans(matrix(as.matrix(colMeans(as.matrix(data))), nrow = 9, ncol = 31)))),A)
}

## --------------------------------------------------------------------------------------------------------------------------------
## Panel B    POT    System → Industries
PBs2i <- matrix(nrow=9,ncol=9)
testall <- read.csv("E:/software/R/code/CoVaRCoes/POTs2ifullCoVaRCoES.csv",head=TRUE)
for (k in 1:3){
  data <- testall[,(279*(k-1)+1):(279*k)]
  result <- foreach(i=1:31,.packages = "MCS")%dopar%{
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
  }
  
  MCSR <- matrix(nrow=9,ncol=31)
  for (j in 1:31){
    r<- result[[j]]
    MCSR[,j] <- r[,2]
  }
  
  for (i in 1:9){
    MCSresult90R[i] <- length(MCSR[i,MCSR[i,]>0.90])
    MCSresult95R[i] <- length(MCSR[i,MCSR[i,]>0.95])
  }
  A=t(rbind(MCSresult95R,MCSresult90R))
  
  PBs2i[,(3*(k-1)+1):(3*k)] <- cbind(t(t(10*rowMeans(matrix(as.matrix(colMeans(as.matrix(data))), nrow = 9, ncol = 31)))),A)
}


##----------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------
##                            Table 4 
##          Panel A Industries → System
PAi2s  

##          Panel A System → Industries 
PAs2i  

##          Panel B Industries → System 
PBi2s 

##          Panel B System → Industries 
PBs2i 
