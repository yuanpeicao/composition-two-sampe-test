#----------------------------------------------------------------
powerFun <- function(x1,x2,propSubset,nRep){
  nRej <- 0
  n <- dim(x1)[1]
  nSubset <- floor(n*propSubset)
  for (i in 1:nRep) {
    index <- sample(n, nSubset, replace = TRUE)
    x1subset <- x1[index,]
    x2subset <- x2[index,]
    pvalue <- pValue(x1subset,x2subset)
    if (pvalue < 0.05) {
      nRej <- nRej + 1
    }
  }
  power <- nRej/nRep
  return(power)
}

#------------------------------------------------------------------------------------
pValue <- function(x1,x2){
  n <- dim(x1)[1]
  p <- ncol(x1)
  SqStandizeDiff <- ((colSums(x1-x2))/n)^2/((diag(var(x1-x2)))*(n-1)/n^2)
  StatX <- max(SqStandizeDiff)
  pvalue <- 1-exp(-1/sqrt(pi)*exp(-(StatX-(2*log(p)-log(log(p))))/2))
  return(pvalue)
}

#------------------------------------------------------------------------------------
pValueTwoSample <- function(x1,x2){
  n1 <- dim(x1)[1]
  n2 <- dim(x2)[1]
  p <- ncol(x1)
  x1mean <- colSums(x1)/n1
  x2mean <- colSums(x2)/n2
  x1Var <- diag(var(x1))*(n1-1)/n1
  x2Var <- diag(var(x2))*(n2-1)/n2
  xStatVar <- (x1Var*n1 + x2Var*n2)/(n1*n2)
  xStat <- max((x1mean/sqrt(xStatVar) - x2mean/sqrt(xStatVar))^2)
  pvalue <- 1-exp(-1/sqrt(pi)*exp(-(xStat-(2*log(p)-log(log(p))))/2))
  return(pvalue)
}
#------------------------------------------------------------------------------------

set.seed(100)

require(ggplot2)
require(reshape)
require(xlsx)
# Load the data
#setwd("dataset")

xt1 <- read.xlsx("please/please-antiTNF1.xlsx", sheetName="Sheet1")
xt2 <- read.xlsx("please/please-antiTNF2.xlsx", sheetName="Sheet1")
xt3 <- read.xlsx("please/please-antiTNF3.xlsx", sheetName="Sheet1")
xt4 <- read.xlsx("please/please-antiTNF4.xlsx", sheetName="Sheet1")

response1 <- xt1[,2]
response2 <- xt2[,2]
response3 <- xt3[,2]
response4 <- xt4[,2]

xt1 <- as.matrix(xt1[,-c(1,2)])
xt2 <- as.matrix(xt2[,-c(1,2)])
xt3 <- as.matrix(xt3[,-c(1,2)])
xt4 <- as.matrix(xt4[,-c(1,2)])

indexResp2 <- which(response2 == 0 & response1 == 0)
indexResp3 <- which(response3 == 0 & response1 == 0)
indexResp4 <- which(response4 == 0 & response1 == 0)

xt12 <- xt1[indexResp2,]
xt13 <- xt1[indexResp3,]
xt14 <- xt1[indexResp4,]
xt2 <- xt2[indexResp2,]
xt3 <- xt3[indexResp3,]
xt4 <- xt4[indexResp4,]

p <- ncol(xt1)

# Filter the data
thresholdt2 <- 1
thresholdt3 <- 2
thresholdt4 <- 2

index12 <- colSums(xt12>0)>thresholdt2 & colSums(xt2>0)>thresholdt2
index13 <- colSums(xt13>0)>thresholdt3 & colSums(xt3>0)>thresholdt3
index14 <- colSums(xt14>0)>thresholdt4 & colSums(xt4>0)>thresholdt4

xt12 <- xt12[,index12]
xt13 <- xt13[,index13]
xt14 <- xt14[,index14]
xt2 <- xt2[,index12]
xt3 <- xt3[,index13]
xt4 <- xt4[,index14]

index <- 1:p
indexT <- index[xTotTB]

# Zero Replacement, replaced by half of the minimum nonzero elements in each treatment
minXt <- min(c(min(xt1[xt1!=0]),min(xt2[xt2!=0]),min(xt3[xt3!=0]),min(xt4[xt4!=0]),min(xt12[xt12!=0]),min(xt13[xt13!=0]),min(xt14[xt14!=0])))

xt12[xt12==0] <- minXt*0.5
xt13[xt13==0] <- minXt*0.5
xt14[xt14==0] <- minXt*0.5
xt2[xt2==0] <- minXt*0.5
xt3[xt3==0] <- minXt*0.5
xt4[xt4==0] <- minXt*0.5

# Covert to the composition
pt2 <- ncol(xt2)
pt3 <- ncol(xt3)
pt4 <- ncol(xt4)
xt12 <- xt12/(rowSums(xt12)%*%matrix(1,1,pt2))
xt13 <- xt13/(rowSums(xt13)%*%matrix(1,1,pt3))
xt14 <- xt14/(rowSums(xt14)%*%matrix(1,1,pt4))
xt2 <- xt2/(rowSums(xt2)%*%matrix(1,1,pt2))
xt3 <- xt3/(rowSums(xt3)%*%matrix(1,1,pt3))
xt4 <- xt4/(rowSums(xt4)%*%matrix(1,1,pt4))


# Logratio and clr
logXt12 <- log(xt12)
logXt13 <- log(xt13)
logXt14 <- log(xt14)
logXt2 <- log(xt2)
logXt3 <- log(xt3)
logXt4 <- log(xt4)

clrXt12 <- logXt12-1/pt2*rowSums(logXt12)%*%matrix(1,1,pt2)
clrXt13 <- logXt13-1/pt3*rowSums(logXt13)%*%matrix(1,1,pt3)
clrXt14 <- logXt14-1/pt4*rowSums(logXt14)%*%matrix(1,1,pt4)
clrXt2 <- logXt2-1/pt2*rowSums(logXt2)%*%matrix(1,1,pt2)
clrXt3 <- logXt3-1/pt3*rowSums(logXt3)%*%matrix(1,1,pt3)
clrXt4 <- logXt4-1/pt4*rowSums(logXt4)%*%matrix(1,1,pt4)

# Paired Test
pXt12 <- pValue(xt12,xt2)
pXt13 <- pValue(xt13,xt3)
pXt14 <- pValue(xt14,xt4)

pLogXt12 <- pValue(logXt12,logXt2)
pLogXt13 <- pValue(logXt13,logXt3)
pLogXt14 <- pValue(logXt14,logXt4)

pclrXt12 <- pValue(clrXt12,clrXt2)
pclrXt13 <- pValue(clrXt13,clrXt3)
pclrXt14 <- pValue(clrXt14,clrXt4)

# empirical power
powerX12 = c()
powerLogX12 = c()
powerClrX12 = c()
for (i in 4:10){
  powerX12 <- c(powerX12,powerFun(xt12,xt2,i*0.1,1000))
  powerLogX12 <- c(powerLogX12,powerFun(logXt12,logXt2,i*0.1,1000))
  powerClrX12 <- c(powerClrX12,powerFun(clrXt12,clrXt2,i*0.1,1000))
}

# plot the curve
grid = c(rep(c(4:10)/10, 3))
power = c(powerClrX12,powerLogX12,powerX12)
Method = c(rep(c("Proposed"), 7),rep(c("Log"), 7),rep(c("Raw"), 7))
EmpiricalPower<-data.frame(Grid = grid, Power = power, Method = Method)
setEPS(width=8, height=6)
postscript(file="C:/Users/yuanpei/Desktop/two sample test/20160916realdata/pleasePower12.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
powerPlot <- qplot(Grid, Power, data = EmpiricalPower, shape = Method, linetype = Method,geom = "path", xlab = "Proportion of sample size", ylab = "Power")
powerPlot + scale_colour_grey() + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(face="bold", 
                                                                                                                                                     size="12"), legend.position="right")+geom_point(size=3)
dev.off()

# empirical power
powerX13 = c()
powerLogX13 = c()
powerClrX13 = c()
for (i in 4:10){
  powerX13 <- c(powerX13,powerFun(xt13,xt3,i*0.1,1000))
  powerLogX13 <- c(powerLogX13,powerFun(logXt13,logXt3,i*0.1,1000))
  powerClrX13 <- c(powerClrX13,powerFun(clrXt13,clrXt3,i*0.1,1000))
}

# plot the curve
grid = c(rep(c(4:10)/10, 3))
power = c(powerClrX13,powerLogX13,powerX13)
Method = c(rep(c("Proposed"), 7),rep(c("Log"), 7),rep(c("Raw"), 7))
EmpiricalPower<-data.frame(Grid = grid, Power = power, Method = Method)
setEPS(width=8, height=6)
postscript(file="C:/Users/yuanpei/Desktop/two sample test/20160916realdata/pleasePower13.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
powerPlot <- qplot(Grid, Power, data = EmpiricalPower, shape = Method, linetype = Method,geom = "path", xlab = "Proportion of sample size", ylab = "Power")
powerPlot + scale_colour_grey() + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(face="bold", 
                                                                                                                                                     size="12"), legend.position="right")+geom_point(size=3)
dev.off()

# empirical power
powerX14 = c()
powerLogX14 = c()
powerClrX14 = c()
for (i in 4:10){
  powerX14 <- c(powerX14,powerFun(xt14,xt4,i*0.1,1000))
  powerLogX14 <- c(powerLogX14,powerFun(logXt14,logXt4,i*0.1,1000))
  powerClrX14 <- c(powerClrX14,powerFun(clrXt14,clrXt4,i*0.1,1000))
}

# plot the curve
grid = c(rep(c(4:10)/10, 3))
power = c(powerClrX14,powerLogX14,powerX14)
Method = c(rep(c("Proposed"), 7),rep(c("Log"), 7),rep(c("Raw"), 7))
EmpiricalPower<-data.frame(Grid = grid, Power = power, Method = Method)
setEPS(width=8, height=6)
postscript(file="C:/Users/yuanpei/Desktop/two sample test/20160916realdata/pleasePower14.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
powerPlot <- qplot(Grid, Power, data = EmpiricalPower, shape = Method, linetype = Method,geom = "path", xlab = "Proportion of sample size", ylab = "Power")
powerPlot + scale_colour_grey() + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(face="bold", 
                                                                                                                                                     size="12"), legend.position="right")+geom_point(size=3)
dev.off()