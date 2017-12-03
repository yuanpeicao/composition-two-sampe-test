#----------------------------------------------------------------
powerFun <- function(x1,x2,propSubset,nRep){
  nRej <- 0
  n1 <- dim(x1)[1]
  n2 <- dim(x2)[1]
  n1Subset <- ceiling(n1*propSubset)
  n2Subset <- ceiling(n2*propSubset)
  for (i in 1:nRep) {
    x1subset <- x1[sample(n1, n1Subset, replace = TRUE),]
    x2subset <- x2[sample(n2, n2Subset, replace = TRUE),]
    pX <- pValue(x1subset,x2subset)
    if (pX < 0.05) {
      nRej <- nRej + 1
    }
  }
  power <- nRej/nRep
  return(power)
}
#------------------------------------------------------------------------------------
pValue <- function(x1,x2){
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
pValueBS <- function(x1,x2,nRep){
  nx1 <- dim(x1)[1]
  nx2 <- dim(x2)[1]
  pValueList <- c()
  for (i in 1:nRep) {
    x1BSSet <- x1[sample(nx1, nx1, replace = TRUE),]
    x2BSSet <- x2[sample(nx2, nx2, replace = TRUE),]
    pValueList <- c(pValueList,pValue(x1BSSet,x2BSSet))
  }
  return(pValueList)
}
#------------------------------------------------------------------------------------
pValueBkTest <- function(x1,x2,nRep){
  n1 <- dim(x1)[1]
  n2 <- dim(x2)[1]
  xTot <- rbind2(x1,x2)
  pvalueList <- c()
  for (i in 1:nRep) {
    index <- sample(n1+n2, n1+n2, replace = FALSE)
    x1BkTSet <- xTot[index[1:n1],]
    x2BkTSet <- xTot[index[(n1+1):(n1+n2)],]
    pvalueList <- c(pvalueList,pValue(x1BkTSet,x2BkTSet))
  }
  return(pvalueList)
}
#------------------------------------------------------------------------------------

require(ggplot2)
require(reshape)

set.seed(350)

# Load the data
w <- as.matrix(read.table("data/bmi_counts.dat"))

# Filter the data
n <- dim(w)[1]

# Add 0.5 on zeros
w[w==0] <- w[w==0]+0.5

# Delete the 58th bacteria (Allisonella)
#w <- w[,-58]
#w <- w[,-36]
p<-dim(w)[2]

# Separate the sample into two groups
quan<-0.25
BMI_and_Counts <- read.table("data/BMI_and_Counts.txt", header=T, quote="\"")
bmi<-BMI_and_Counts[,2]

index_1<-which(bmi>=25)
index_2<-which(bmi<25)

wx <- w[index_1,]
wy <- w[index_2,]

nx <- dim(wx)[1]
ny <- dim(wy)[1]

x <- wx/(rowSums(wx)%*%matrix(1,1,p))
y <- wy/(rowSums(wy)%*%matrix(1,1,p))

log_x <- log(x)
log_y <- log(y)

clog_TX <- log_x-1/p*rowSums(log_x)%*%matrix(1,1,p)
clog_TY <- log_y-1/p*rowSums(log_y)%*%matrix(1,1,p)

# p-values
p_x <- pValue(x,y)
p_logx <- pValue(log_x,log_y)
p_clrx <- pValue(clog_TX,clog_TY)

# empirical power, calculate the power on the grid 0.2:0.05:0.1
powerX <- c()
powerLogX <- c()
powerClrX <- c()
for (i in 4:20){
  powerX <- c(powerX,powerFun(x,y,i*0.05,5000))
  powerLogX <- c(powerLogX,powerFun(log_x,log_y,i*0.05,5000))
  powerClrX <- c(powerClrX,powerFun(clog_TX,clog_TY,i*0.05,5000))
}

# plot the curve
grid = c(rep(c(4:20)/20, 3)) # set the grid 0.2 : 0.05 : 1
power = c(powerClrX,powerLogX,powerX)
Method = c(rep(c("Proposed"), 17),rep(c("Log"), 17),rep(c("Raw"), 17))
EmpiricalPower<-data.frame(Grid = grid, Power = power, Method = Method)
setEPS(width=8, height=6)
postscript(file="Plot/comboPower.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
powerPlot <- qplot(Grid, Power, data = EmpiricalPower, shape = Method, linetype = Method,geom = "path", xlab = "Proportion of sample size", ylab = "Power")
powerPlot + scale_colour_grey() + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(face="bold", 
                                                             size="12"), legend.position="right")+geom_point(size=3)
dev.off()

# histogram for bootstrap resampling
pValueBSList <- pValueBS(clog_TX,clog_TY,1000)
setEPS(width=10, height=6)
postscript(file="Plot/comboBS.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
hist(pValueBSList, xlim=c(0,0.4), ylim=c(0,1000),breaks = 20, xlab = "p-value", ylab = "Frequency", main = NULL, cex = 2)
dev.off()

# histogram for backtesting bootstrap resampling
pValueBkTestList <- pValueBkTest(clog_TX,clog_TY,1000)
setEPS(width=10, height=6)
postscript(file="C:/Users/yuanpei/Desktop/two sample test/20160916realdata/comboBkTest.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
hist(pValueBkTestList, xlim=c(0,1.0), ylim=c(0,200), breaks = 20, xlab = "p-value", ylab = "Frequency", main = NULL, cex = 2)
dev.off()