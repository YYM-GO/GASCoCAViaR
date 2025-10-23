
rm(list = ls())

library(superheat)
setwd("E:/software/R/code/CoVaRCoes")
source("E:/software/R/code/CoVaRCoes/trafficlight.R")

##-----------------------------------------------------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------------------------------------------
##        (a) Model comparison: Industry → System  Gaussian innovations
# (1) left top
data <- read.csv("i2sOOSGauSAVAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
          "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: Industry to System (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (2) right top
data <- read.csv("i2sOOSGauASAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
          "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: Industry to System (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (3) left middle
data <- read.csv("i2sOOSGauSAVN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
          "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: Industry to System (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (4) right middle
data <- read.csv("i2sOOSGauASN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
          "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: Industry to System (FHN)"
trafficlight(data,name,plottitle )   #(result)


# (5) left bottom
data <- read.csv("i2sOOSGauSAV0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
          "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: Industry to System (FH0)"
trafficlight(data,name,plottitle )   #(result)

# (6) right bottom
data <- read.csv("i2sOOSGauAS0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
          "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: Industry to System (FH0)"
trafficlight(data,name,plottitle )   #(result)


##----------------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------
##        (b) Model comparison: System → Industry  Gaussian innovations

# (1) left top
data <- read.csv("s2iOOSGauSAVAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
          "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: System to Industry (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (2) right top
data <- read.csv("s2iOOSGauASAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
          "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: System to Industry (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (3) left middle
data <- read.csv("s2iOOSGauSAVN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
          "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: System to Industry (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (4) right middle
data <- read.csv("s2iOOSGauASN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
          "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: System to Industry (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (5) left bottom
data <- read.csv("s2iOOSGauSAV0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
          "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: System to Industry (FH0)"
trafficlight(data,name,plottitle )   #(result)

# (6) right bottom
data <- read.csv("s2iOOSGauAS0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
          "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: System to Industry (FH0)"
trafficlight(data,name,plottitle )   #(result)







##----------------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------
##            (c) Model comparison: Industry → System     POT innovations
# (1) left top
data <- read.csv("i2sOOSPOTSAVAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
          "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: Industry to System (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (2) right top
data <- read.csv("i2sOOSPOTASAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
          "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: Industry to System (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (3) left middle
data <- read.csv("i2sOOSPOTSAVN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
          "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: Industry to System (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (4) right middle
data <- read.csv("i2sOOSPOTASN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
          "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: Industry to System (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (5) left bottom
data <- read.csv("i2sOOSPOTSAV0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
          "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: Industry to System (FH0)"
trafficlight(data,name,plottitle )   #(result)

# (6) right bottom
data <- read.csv("i2sOOSPOTAS0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
          "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: Industry to System (FH0)"
trafficlight(data,name,plottitle )   #(result)

##----------------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------
##         (d) Model comparison: System → Industry    POT innovations
# (1) left top
data <- read.csv("s2iOOSPOTSAVAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
          "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: System to Industry (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (2) right top
data <- read.csv("s2iOOSPOTASAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
          "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: System to Industry (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (3) left middle
data <- read.csv("s2iOOSPOTSAVN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
          "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: System to Industry (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (4) right middle
data <- read.csv("s2iOOSPOTASN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
          "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: System to Industry (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (5) left bottom
data <- read.csv("s2iOOSPOTSAV0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
          "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: System to Industry (FH0)"
trafficlight(data,name,plottitle )   #(result)

# (6) right bottom
data <- read.csv("s2iOOSPOTAS0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
          "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: System to Industry (FH0)"
trafficlight(data,name,plottitle )   #(result)