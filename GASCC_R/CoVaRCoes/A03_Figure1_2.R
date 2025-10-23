    
rm(list = ls())

library(superheat)
setwd("E:/software/R/code/CoVaRCoes")
source("E:/software/R/code/CoVaRCoes/trafficlight.R")

##-----------------------------------------------------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------------------------------------------
##        (a) Model comparison: Industry → System  Gaussian innovations
# (1) left top
data <- read.csv("i2sfullGauSAVAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
     "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: Industry to System (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (2) right top
data <- read.csv("i2sfullGauASAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
                    "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: Industry to System (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (3) left middle
data <- read.csv("i2sfullGauSAVN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
                    "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: Industry to System (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (4) right middle
data <- read.csv("i2sfullGauASN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
                    "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: Industry to System (FHN)"
trafficlight(data,name,plottitle )   #(result)


# (5) left bottom
data <- read.csv("i2sfullGauSAV0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
                    "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: Industry to System (FH0)"
trafficlight(data,name,plottitle )   #(result)

# (6) right bottom
data <- read.csv("i2sfullGauAS0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
                    "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: Industry to System (FH0)"
trafficlight(data,name,plottitle )   #(result)


##----------------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------
##        (b) Model comparison: System → Industry  Gaussian innovations

# (1) left top
data <- read.csv("s2ifullGauSAVAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
                    "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: System to Industry (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (2) right top
data <- read.csv("s2ifullGauASAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
                    "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: System to Industry (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (3) left middle
data <- read.csv("s2ifullGauSAVN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
                    "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: System to Industry (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (4) right middle
data <- read.csv("s2ifullGauASN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
                    "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: System to Industry (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (5) left bottom
data <- read.csv("s2ifullGauSAV0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
                    "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: System to Industry (FH0)"
trafficlight(data,name,plottitle )   #(result)

# (6) right bottom
data <- read.csv("s2ifullGauAS0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
                    "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: System to Industry (FH0)"
trafficlight(data,name,plottitle )   #(result)







##----------------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------
##            (c) Model comparison: Industry → System     POT innovations
# (1) left top
data <- read.csv("i2sfullPOTSAVAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
                    "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: Industry to System (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (2) right top
data <- read.csv("i2sfullPOTASAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
                    "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: Industry to System (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (3) left middle
data <- read.csv("i2sfullPOTSAVN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
                    "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: Industry to System (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (4) right middle
data <- read.csv("i2sfullPOTASN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
                    "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: Industry to System (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (5) left bottom
data <- read.csv("i2sfullPOTSAV0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
                    "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: Industry to System (FH0)"
trafficlight(data,name,plottitle )   #(result)

# (6) right bottom
data <- read.csv("i2sfullPOTAS0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
                    "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: Industry to System (FH0)"
trafficlight(data,name,plottitle )   #(result)

##----------------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------
##         (d) Model comparison: System → Industry    POT innovations
# (1) left top
data <- read.csv("s2ifullPOTSAVAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
                    "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: System to Industry (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (2) right top
data <- read.csv("s2ifullPOTASAL.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
                    "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: System to Industry (FHAL)"
trafficlight(data,name,plottitle )   #(result)

# (3) left middle
data <- read.csv("s2ifullPOTSAVN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
                    "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: System to Industry (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (4) right middle
data <- read.csv("s2ifullPOTASN.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
                    "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: System to Industry (FHN)"
trafficlight(data,name,plottitle )   #(result)

# (5) left bottom
data <- read.csv("s2ifullPOTSAV0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoSAV:CoSAVc","GAS-CoSAV:GAS-CoAS","GAS-CoSAV:CoASc",
                    "GAS-CoSAV:Gaussian","GAS-CoSAV:Student's-t","GAS-CoSAV:Clayton","GAS-CoSAV:RGumbel","GAS-CoSAV:DCC")
plottitle <- " GAS-CoSAV: System to Industry (FH0)"
trafficlight(data,name,plottitle )   #(result)

# (6) right bottom
data <- read.csv("s2ifullPOTAS0.csv",header=F)
data <- as.data.frame(data)
name <- c("GAS-CoAS:CoASc","GAS-CoAS:GAS-CoSAV","GAS-CoAS:CoSAVc",
                    "GAS-CoAS:Gaussian","GAS-CoAS:Student's-t","GAS-CoAS:Clayton","GAS-CoAS:RGumbel","GAS-CoAS:DCC")
plottitle <- " GAS-CoAS: System to Industry (FH0)"
trafficlight(data,name,plottitle )   #(result)