
rm(list = ls())
setwd("E:/software/R/code/CoVaRCoes")
library("igraph") 


nodes <- read.csv("vertex.csv", header=T, as.is=T)

##  low risk  (left)
links <- read.csv("sysriskgood.csv", header=T, as.is=T)
head(nodes)

head(links)

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
colrs1 <- rainbow(3, alpha=.5)
colrs10 <- rainbow(3, alpha=1)
colrs2 <- c("#DF65B0", "#252525")
V(net)$color <- colrs1[V(net)$media.type]
V(net)$frame.color <- colrs10[V(net)$media.type]
E(net)$edge.color <- colrs2[E(net)$type]

E(net)$arrow.mode <- 0
edge.start <- ends(net, es=E(net), names=F)[,1]
edge.col <- E(net)$color2[edge.start]
position<-c(0,0,0,0,-pi/5,-pi/4,-pi/3,-pi/2,-pi/2,-pi/1.6,-pi/1.5,pi/0.8,pi/0.9,pi,pi,pi,pi,pi,pi,pi/1.1,pi/1.2,pi/1.5,pi/2,pi/2,pi/2,pi/3,pi/4,pi/5,pi/6,0,0)

l <- layout_in_circle(net)
plot(net, layout=l,edge.color=E(net)$edge.color, edge.curved=0,vertex.size=15,vertex.label=V(net)$media,
     vertex.label.dist=2,vertex.label.degree=position,edge.width=1.5)
legend(x="bottomright", c("Upstream","Midstream", "Downstream"), pch=21,col=colrs1, pt.bg=colrs1, 
       pt.cex=2, cex=1, bty="n", ncol=1,xpd = TRUE)



##     high risk  (right)

links <- read.csv("sysriskbad.csv", header=T, as.is=T)
head(nodes)

head(links)
library('igraph')
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
colrs1 <- rainbow(3, alpha=.5)
colrs10 <- rainbow(3, alpha=1)
colrs2 <- c("#DF65B0", "#252525")
V(net)$color <- colrs1[V(net)$media.type]
V(net)$frame.color <- colrs10[V(net)$media.type]
E(net)$edge.color <- colrs2[E(net)$type]

E(net)$arrow.mode <- 0
edge.start <- ends(net, es=E(net), names=F)[,1]
edge.col <- E(net)$color2[edge.start]
position<-c(0,0,0,0,-pi/5,-pi/4,-pi/3,-pi/2,-pi/2,-pi/1.6,-pi/1.5,pi/0.8,pi/0.9,pi,pi,pi,pi,pi,pi,pi/1.1,pi/1.2,pi/1.5,pi/2,pi/2,pi/2,pi/3,pi/4,pi/5,pi/6,0,0)

l <- layout_in_circle(net)


plot(net, layout=l,edge.color=E(net)$edge.color, edge.curved=0,vertex.size=15,vertex.label=V(net)$media,
     vertex.label.dist=2,vertex.label.degree=position,edge.width=1.5)

legend(x="bottomright", c("Upstream","Midstream", "Downstream"), pch=21,col=colrs1, pt.bg=colrs1, 
       pt.cex=2, cex=1, bty="n", ncol=1,xpd = TRUE)
#title("Network visualization on February 3, 2020")
