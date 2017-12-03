require(ggplot2)
require(reshape)
require(xlsx)
# Load the data
#setwd("dataset")

xt1 <- read.xlsx("please/please-antiTNF1.xlsx", sheetName="Sheet1")
xt2 <- read.xlsx("please/please-antiTNF2.xlsx", sheetName="Sheet1")
xt3 <- read.xlsx("please/please-antiTNF3.xlsx", sheetName="Sheet1")
xt4 <- read.xlsx("please/please-antiTNF4.xlsx", sheetName="Sheet1")

xt1 <- xt1[,-1]
xt2 <- xt2[,-1]
xt3 <- xt3[,-1]
xt4 <- xt4[,-1]

# deleted selected genera based on clr
#xt1 <- xt1[,-which(names(xt1) %in% c("Neisseria","Phascolarctobacterium","Candidatus???Zinderia","Fusobacterium","Subdoligranulum","Methanobrevibacter","Anaerostipes","Shigella","Staphylococcus","Blautia"))]
#xt2 <- xt2[,-which(names(xt2) %in% c("Neisseria","Phascolarctobacterium","Candidatus???Zinderia","Fusobacterium","Subdoligranulum","Methanobrevibacter","Anaerostipes","Shigella","Staphylococcus","Blautia"))]
#xt3 <- xt3[,-which(names(xt3) %in% c("Neisseria","Phascolarctobacterium","Candidatus???Zinderia","Fusobacterium","Subdoligranulum","Methanobrevibacter","Anaerostipes","Shigella","Staphylococcus","Blautia"))]
#xt4 <- xt4[,-which(names(xt4) %in% c("Neisseria","Phascolarctobacterium","Candidatus???Zinderia","Fusobacterium","Subdoligranulum","Methanobrevibacter","Anaerostipes","Shigella","Staphylococcus","Blautia"))]

# deleted selected genera based on LogX
#xt1 <- xt1[,-which(names(xt1) %in% c("Phascolarctobacterium","Methanobrevibacter","Neisseria","Candidatus???Zinderia","Fusobacterium","Subdoligranulum","Staphylococcus","Shigella","Anaerofustis","Eikenella"))]
#xt2 <- xt2[,-which(names(xt2) %in% c("Phascolarctobacterium","Methanobrevibacter","Neisseria","Candidatus???Zinderia","Fusobacterium","Subdoligranulum","Staphylococcus","Shigella","Anaerofustis","Eikenella"))]
#xt3 <- xt3[,-which(names(xt3) %in% c("Phascolarctobacterium","Methanobrevibacter","Neisseria","Candidatus???Zinderia","Fusobacterium","Subdoligranulum","Staphylococcus","Shigella","Anaerofustis","Eikenella"))]
#xt4 <- xt4[,-which(names(xt4) %in% c("Phascolarctobacterium","Methanobrevibacter","Neisseria","Candidatus???Zinderia","Fusobacterium","Subdoligranulum","Staphylococcus","Shigella","Anaerofustis","Eikenella"))]

xt1 <- as.matrix(xt1)
xt2 <- as.matrix(xt2)
xt3 <- as.matrix(xt3)
xt4 <- as.matrix(xt4)

xd1 <- read.xlsx("please/please-diet1.xlsx", sheetName="Sheet1")
xd2 <- read.xlsx("please/please-diet2.xlsx", sheetName="Sheet1")
xd3 <- read.xlsx("please/please-diet3.xlsx", sheetName="Sheet1")
xd4 <- read.xlsx("please/please-diet4.xlsx", sheetName="Sheet1")

xd1 <- xd1[,-1]
xd2 <- xd2[,-1]
xd3 <- xd3[,-1]
xd4 <- xd4[,-1]

# deleted selected genera based on clr
#xd1 <- xd1[,-which(names(xd1) %in% c("Abiotrophia","Anaerostipes","Haemophilus"))]
#xd2 <- xd2[,-which(names(xd2) %in% c("Abiotrophia","Anaerostipes","Haemophilus"))]
#xd3 <- xd3[,-which(names(xd3) %in% c("Abiotrophia","Anaerostipes","Haemophilus"))]
#xd4 <- xd4[,-which(names(xd4) %in% c("Abiotrophia","Anaerostipes","Haemophilus"))]

# deleted selected genera based on LogX
#xd1 <- xd1[,-which(names(xd1) %in% c("Anaerostipes","Abiotrophia","Proteus","Shigella","Haemophilus"))]
#xd2 <- xd2[,-which(names(xd2) %in% c("Anaerostipes","Abiotrophia","Proteus","Shigella","Haemophilus"))]
#xd3 <- xd3[,-which(names(xd3) %in% c("Anaerostipes","Abiotrophia","Proteus","Shigella","Haemophilus"))]
#xd4 <- xd4[,-which(names(xd4) %in% c("Anaerostipes","Abiotrophia","Proteus","Shigella","Haemophilus"))]

xd1 <- as.matrix(xd1)
xd2 <- as.matrix(xd2)
xd3 <- as.matrix(xd3)
xd4 <- as.matrix(xd4)

p <- ncol(xt1)

# Filter the data
threshold = 4
xt1b <- colSums(xt1>0)>threshold
xt2b <- colSums(xt2>0)>threshold
xt3b <- colSums(xt3>0)>threshold
xt4b <- colSums(xt4>0)>threshold
xd1b <- colSums(xd1>0)>threshold
xd2b <- colSums(xd2>0)>threshold
xd3b <- colSums(xd3>0)>threshold
xd4b <- colSums(xd4>0)>threshold

xTotTB <- xt1b & xt2b & xt3b & xt4b
xTotDB <- xd1b & xd2b & xd3b & xd4b

xt1 <- xt1[,xTotTB]
xt2 <- xt2[,xTotTB]
xt3 <- xt3[,xTotTB]
xt4 <- xt4[,xTotTB]

xd1 <- xd1[,xTotDB]
xd2 <- xd2[,xTotDB]
xd3 <- xd3[,xTotDB]
xd4 <- xd4[,xTotDB]

index <- 1:p
indexT <- index[xTotTB]
indexD <- index[xTotDB]

# Zero Replacement
minXt1 <- min(xt1[xt1!=0])
minXt2 <- min(xt2[xt2!=0])
minXt3 <- min(xt3[xt3!=0])
minXt4 <- min(xt4[xt4!=0])

minXd1 <- min(xd1[xd1!=0])
minXd2 <- min(xd2[xd2!=0])
minXd3 <- min(xd3[xd3!=0])
minXd4 <- min(xd4[xd4!=0])

xt1[xt1==0] <- minXt1*0.5
xt2[xt2==0] <- minXt2*0.5
xt3[xt3==0] <- minXt3*0.5
xt4[xt4==0] <- minXt4*0.5

xd1[xd1==0] <- minXd1*0.5
xd2[xd2==0] <- minXd2*0.5
xd3[xd3==0] <- minXd3*0.5
xd4[xd4==0] <- minXd4*0.5

# Covert to the composition
pt <- ncol(xt1)
xt1 <- xt1/(rowSums(xt1)%*%matrix(1,1,pt))
xt2 <- xt2/(rowSums(xt2)%*%matrix(1,1,pt))
xt3 <- xt3/(rowSums(xt3)%*%matrix(1,1,pt))
xt4 <- xt4/(rowSums(xt4)%*%matrix(1,1,pt))

pd <- ncol(xd1)
xd1 <- xd1/(rowSums(xd1)%*%matrix(1,1,pd))
xd2 <- xd2/(rowSums(xd2)%*%matrix(1,1,pd))
xd3 <- xd3/(rowSums(xd3)%*%matrix(1,1,pd))
xd4 <- xd4/(rowSums(xd4)%*%matrix(1,1,pd))

# Logratio and clr
logXt1 <- log(xt1)
logXt2 <- log(xt2)
logXt3 <- log(xt3)
logXt4 <- log(xt4)

logXd1 <- log(xd1)
logXd2 <- log(xd2)
logXd3 <- log(xd3)
logXd4 <- log(xd4)

clrXt1 <- logXt1-1/pt*rowSums(logXt1)%*%matrix(1,1,pt)
clrXt2 <- logXt2-1/pt*rowSums(logXt2)%*%matrix(1,1,pt)
clrXt3 <- logXt3-1/pt*rowSums(logXt3)%*%matrix(1,1,pt)
clrXt4 <- logXt4-1/pt*rowSums(logXt4)%*%matrix(1,1,pt)

clrXd1 <- logXd1-1/pd*rowSums(logXd1)%*%matrix(1,1,pd)
clrXd2 <- logXd2-1/pd*rowSums(logXd2)%*%matrix(1,1,pd)
clrXd3 <- logXd3-1/pd*rowSums(logXd3)%*%matrix(1,1,pd)
clrXd4 <- logXd4-1/pd*rowSums(logXd4)%*%matrix(1,1,pd)

# Squared standardized difference
nt <- nrow(xt1)
SqStandizeDiffXt12 <- ((colSums(xt1-xt2))/nt)^2/((diag(var(xt1-xt2)))*(nt-1)/nt^2)
SqStandizeDiffXt23 <- ((colSums(xt2-xt3))/nt)^2/((diag(var(xt2-xt3)))*(nt-1)/nt^2)
SqStandizeDiffXt34 <- ((colSums(xt3-xt4))/nt)^2/((diag(var(xt3-xt4)))*(nt-1)/nt^2)

SqStandizeDiffLogXt12 <- ((colSums(logXt1-logXt2))/nt)^2/((diag(var(logXt1-logXt2)))*(nt-1)/nt^2)
SqStandizeDiffLogXt23 <- ((colSums(logXt2-logXt3))/nt)^2/((diag(var(logXt2-logXt3)))*(nt-1)/nt^2)
SqStandizeDiffLogXt34 <- ((colSums(logXt3-logXt4))/nt)^2/((diag(var(logXt3-logXt4)))*(nt-1)/nt^2)

SqStandizeDiffclrXt12 <- ((colSums(clrXt1-clrXt2))/nt)^2/((diag(var(clrXt1-clrXt2)))*(nt-1)/nt^2)
SqStandizeDiffclrXt23 <- ((colSums(clrXt2-clrXt3))/nt)^2/((diag(var(clrXt2-clrXt3)))*(nt-1)/nt^2)
SqStandizeDiffclrXt34 <- ((colSums(clrXt3-clrXt4))/nt)^2/((diag(var(clrXt3-clrXt4)))*(nt-1)/nt^2)

nd <- nrow(xd1)
SqStandizeDiffXd12 <- ((colSums(xd1-xd2))/nd)^2/((diag(var(xd1-xd2)))*(nd-1)/nd^2)
SqStandizeDiffXd23 <- ((colSums(xd2-xd3))/nd)^2/((diag(var(xd2-xd3)))*(nd-1)/nd^2)
SqStandizeDiffXd34 <- ((colSums(xd3-xd4))/nd)^2/((diag(var(xd3-xd4)))*(nd-1)/nd^2)

SqStandizeDiffLogXd12 <- ((colSums(logXd1-logXd2))/nd)^2/((diag(var(logXd1-logXd2)))*(nd-1)/nd^2)
SqStandizeDiffLogXd23 <- ((colSums(logXd2-logXd3))/nd)^2/((diag(var(logXd2-logXd3)))*(nd-1)/nd^2)
SqStandizeDiffLogXd34 <- ((colSums(logXd3-logXd4))/nd)^2/((diag(var(logXd3-logXd4)))*(nd-1)/nd^2)

SqStandizeDiffclrXd12 <- ((colSums(clrXd1-clrXd2))/nd)^2/((diag(var(clrXd1-clrXd2)))*(nd-1)/nd^2)
SqStandizeDiffclrXd23 <- ((colSums(clrXd2-clrXd3))/nd)^2/((diag(var(clrXd2-clrXd3)))*(nd-1)/nd^2)
SqStandizeDiffclrXd34 <- ((colSums(clrXd3-clrXd4))/nd)^2/((diag(var(clrXd3-clrXd4)))*(nd-1)/nd^2)

# p-values
StatXt12 <- max(SqStandizeDiffXt12)
StatXt23 <- max(SqStandizeDiffXt23)
StatXt34 <- max(SqStandizeDiffXt34)

StatLogXt12 <- max(SqStandizeDiffLogXt12)
StatLogXt23 <- max(SqStandizeDiffLogXt23)
StatLogXt34 <- max(SqStandizeDiffLogXt34)

StatclrXt12 <- max(SqStandizeDiffclrXt12)
StatclrXt23 <- max(SqStandizeDiffclrXt23)
StatclrXt34 <- max(SqStandizeDiffclrXt34)

StatXd12 <- max(SqStandizeDiffXd12)
StatXd23 <- max(SqStandizeDiffXd23)
StatXd34 <- max(SqStandizeDiffXd34)

StatLogXd12 <- max(SqStandizeDiffLogXd12)
StatLogXd23 <- max(SqStandizeDiffLogXd23)
StatLogXd34 <- max(SqStandizeDiffLogXd34)

StatclrXd12 <- max(SqStandizeDiffclrXd12)
StatclrXd23 <- max(SqStandizeDiffclrXd23)
StatclrXd34 <- max(SqStandizeDiffclrXd34)

pXt12 <- 1-exp(-1/sqrt(pi)*exp(-(StatXt12-(2*log(pt)-log(log(pt))))/2))
pXt23 <- 1-exp(-1/sqrt(pi)*exp(-(StatXt23-(2*log(pt)-log(log(pt))))/2))
pXt34 <- 1-exp(-1/sqrt(pi)*exp(-(StatXt34-(2*log(pt)-log(log(pt))))/2))

pLogXt12 <- 1-exp(-1/sqrt(pi)*exp(-(StatLogXt12-(2*log(pt)-log(log(pt))))/2))
pLogXt23 <- 1-exp(-1/sqrt(pi)*exp(-(StatLogXt23-(2*log(pt)-log(log(pt))))/2))
pLogXt34 <- 1-exp(-1/sqrt(pi)*exp(-(StatLogXt34-(2*log(pt)-log(log(pt))))/2))

pclrXt12 <- 1-exp(-1/sqrt(pi)*exp(-(StatclrXt12-(2*log(pt)-log(log(pt))))/2))
pclrXt23 <- 1-exp(-1/sqrt(pi)*exp(-(StatclrXt23-(2*log(pt)-log(log(pt))))/2))
pclrXt34 <- 1-exp(-1/sqrt(pi)*exp(-(StatclrXt34-(2*log(pt)-log(log(pt))))/2))

pXd12 <- 1-exp(-1/sqrt(pi)*exp(-(StatXd12-(2*log(pd)-log(log(pd))))/2))
pXd23 <- 1-exp(-1/sqrt(pi)*exp(-(StatXd23-(2*log(pd)-log(log(pd))))/2))
pXd34 <- 1-exp(-1/sqrt(pi)*exp(-(StatXd34-(2*log(pd)-log(log(pd))))/2))

pLogXd12 <- 1-exp(-1/sqrt(pi)*exp(-(StatLogXd12-(2*log(pd)-log(log(pd))))/2))
pLogXd23 <- 1-exp(-1/sqrt(pi)*exp(-(StatLogXd23-(2*log(pd)-log(log(pd))))/2))
pLogXd34 <- 1-exp(-1/sqrt(pi)*exp(-(StatLogXd34-(2*log(pd)-log(log(pd))))/2))

pclrXd12 <- 1-exp(-1/sqrt(pi)*exp(-(StatclrXd12-(2*log(pd)-log(log(pd))))/2))
pclrXd23 <- 1-exp(-1/sqrt(pi)*exp(-(StatclrXd23-(2*log(pd)-log(log(pd))))/2))
pclrXd34 <- 1-exp(-1/sqrt(pi)*exp(-(StatclrXd34-(2*log(pd)-log(log(pd))))/2))

# Load genus and phyla
Genera <- as.matrix(read.table("please/please-genus.dat"))
GeneraT <- Genera[indexT]
GeneraD <- Genera[indexD]


Phylum <- as.matrix(read.table("please/please-phylum.dat"))
PhylumT <- Phylum[indexT]
PhylumD <- Phylum[indexD]

# Selected genus
GenusXt12 <- indexT[which(SqStandizeDiffXt12 == StatXt12)]
GenusXt23 <- indexT[which(SqStandizeDiffXt23 == StatXt23)]
GenusXt34 <- indexT[which(SqStandizeDiffXt34 == StatXt34)]

GenusLogXt12 <- indexT[which(SqStandizeDiffLogXt12 == StatLogXt12)]
GenusLogXt23 <- indexT[which(SqStandizeDiffLogXt23 == StatLogXt23)]
GenusLogXt34 <- indexT[which(SqStandizeDiffLogXt34 == StatLogXt34)]

GenusclrXt12 <- indexT[which(SqStandizeDiffclrXt12 == StatclrXt12)]
GenusclrXt23 <- indexT[which(SqStandizeDiffclrXt23 == StatclrXt23)]
GenusclrXt34 <- indexT[which(SqStandizeDiffclrXt34 == StatclrXt34)]

GenusXd12 <- indexD[which(SqStandizeDiffXd12 == StatXd12)]
GenusXd23 <- indexD[which(SqStandizeDiffXd23 == StatXd23)]
GenusXd34 <- indexD[which(SqStandizeDiffXd34 == StatXd34)]

GenusLogXd12 <- indexD[which(SqStandizeDiffLogXd12 == StatLogXd12)]
GenusLogXd23 <- indexD[which(SqStandizeDiffLogXd23 == StatLogXd23)]
GenusLogXd34 <- indexD[which(SqStandizeDiffLogXd34 == StatLogXd34)]

GenusclrXd12 <- indexD[which(SqStandizeDiffclrXd12 == StatclrXd12)]
GenusclrXd23 <- indexD[which(SqStandizeDiffclrXd23 == StatclrXd23)]
GenusclrXd34 <- indexD[which(SqStandizeDiffclrXd34 == StatclrXd34)]

GNameXt12 <- GeneraT[which(SqStandizeDiffXt12 == StatXt12)]
GNameXt23 <- GeneraT[which(SqStandizeDiffXt23 == StatXt23)]
GNameXt34 <- GeneraT[which(SqStandizeDiffXt34 == StatXt34)]
GNameLogXt12 <- GeneraT[which(SqStandizeDiffLogXt12 == StatLogXt12)]
GNameLogXt23 <- GeneraT[which(SqStandizeDiffLogXt23 == StatLogXt23)]
GNameLogXt34 <- GeneraT[which(SqStandizeDiffLogXt34 == StatLogXt34)]
GNameclrXt12 <- GeneraT[which(SqStandizeDiffclrXt12 == StatclrXt12)]
GNameclrXt23 <- GeneraT[which(SqStandizeDiffclrXt23 == StatclrXt23)]
GNameclrXt34 <- GeneraT[which(SqStandizeDiffclrXt34 == StatclrXt34)]

GNameXd12 <- GeneraD[which(SqStandizeDiffXd12 == StatXd12)]
GNameXd23 <- GeneraD[which(SqStandizeDiffXd23 == StatXd23)]
GNameXd34 <- GeneraD[which(SqStandizeDiffXd34 == StatXd34)]
GNameLogXd12 <- GeneraD[which(SqStandizeDiffLogXd12 == StatLogXd12)]
GNameLogXd23 <- GeneraD[which(SqStandizeDiffLogXd23 == StatLogXd23)]
GNameLogXd34 <- GeneraD[which(SqStandizeDiffLogXd34 == StatLogXd34)]
GNameclrXd12 <- GeneraD[which(SqStandizeDiffclrXd12 == StatclrXd12)]
GNameclrXd23 <- GeneraD[which(SqStandizeDiffclrXd23 == StatclrXd23)]
GNameclrXd34 <- GeneraD[which(SqStandizeDiffclrXd34 == StatclrXd34)]

## Barplot
# antiTNF
# 1 vs 2 (clrX)
clrxT12 <- data.frame(Genera = GeneraT, clrX = SqStandizeDiffclrXt12,Phylum = PhylumT)
clrxT12 <- clrxT12[order(clrxT12$clrX,decreasing = TRUE),]
clrxT12 <- clrxT12[1:10,]
setEPS(width=8, height=4)
postscript(file="please-tnf12-clr.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
p <- ggplot(clrxT12, aes(x=reorder(Genera, clrX),y=clrX,fill=Phylum))
p <- p + geom_bar(stat="identity", position = "identity", width=0.5)
p <- p + labs(x="", y="Squared standardized difference on Y")
p <- p + scale_y_continuous(breaks=seq(0,25,by=2))+ theme(panel.grid.major=element_line(colour=NA), text = element_text(size=15), axis.text.x = element_text(size  = 15, angle = 90, hjust = 1, vjust = 1))
p <- p + coord_flip()
p
dev.off()

# 1 vs 2 (logX)
LogxT12 <- data.frame(Genera = GeneraT, LogX = SqStandizeDiffLogXt12,Phylum = PhylumT)
LogxT12 <- LogxT12[order(LogxT12$LogX,decreasing = TRUE),]
LogxT12 <- LogxT12[1:10,]
setEPS(width=8, height=4)
postscript(file="please-tnf12-Log.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
p <- ggplot(LogxT12, aes(x=reorder(Genera, LogX),y=LogX,fill=Phylum))
p <- p + geom_bar(stat="identity", position = "identity", width=0.5)
p <- p + labs(x="", y="Squared standardized difference on LogX")
p <- p + scale_y_continuous(breaks=seq(0,6000,by=600))+ theme(panel.grid.major=element_line(colour=NA), text = element_text(size=15), axis.text.x = element_text(size  = 15, angle = 90, hjust = 1, vjust = 1))
p <- p + coord_flip()
p
dev.off()

# Diet
# 1 vs 2 (clrX)
clrxD12 <- data.frame(Genera = GeneraD, clrX = SqStandizeDiffclrXd12,Phylum = PhylumD)
clrxD12 <- clrxD12[order(clrxD12$clrX,decreasing = TRUE),]
clrxD12 <- clrxD12[1:10,]
setEPS(width=8, height=4)
postscript(file="please-diet12-clr.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
p <- ggplot(clrxD12, aes(x=reorder(Genera, clrX),y=clrX,fill=Phylum))
p <- p + geom_bar(stat="identity", position = "identity", width=0.5)
p <- p + labs(x="", y="Squared standardized difference on Y")
p <- p + scale_y_continuous(breaks=seq(0,25,by=2))+ theme(panel.grid.major=element_line(colour=NA), text = element_text(size=15), axis.text.x = element_text(size  = 15, angle = 90, hjust = 1, vjust = 1))
p <- p + coord_flip()
p
dev.off()

# 1 vs 2 (LogX)
LogxD12 <- data.frame(Genera = GeneraD, LogX = SqStandizeDiffLogXd12,Phylum = PhylumD)
LogxD12 <- LogxD12[order(LogxD12$LogX,decreasing = TRUE),]
LogxD12 <- LogxD12[1:10,]
setEPS(width=8, height=4)
postscript(file="please-diet12-Log.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
p <- ggplot(LogxD12, aes(x=reorder(Genera, LogX),y=LogX,fill=Phylum))
p <- p + geom_bar(stat="identity", position = "identity", width=0.5)
p <- p + labs(x="", y="Squared standardized difference on LogX")
p <- p + scale_y_continuous(breaks=seq(0,500,by=50))+ theme(panel.grid.major=element_line(colour=NA), text = element_text(size=15), axis.text.x = element_text(size  = 15, angle = 90, hjust = 1, vjust = 1))
p <- p + coord_flip()
p
dev.off()