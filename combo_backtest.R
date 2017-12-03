#----------------------------------------------------------------
powerFun <- function(x1,x2,propSubset,nRep){
  nRej <- 0
  nx1 <- dim(x1)[1]
  nx2 <- dim(x2)[1]
  nx1Subset <- floor(nx1*propSubset)
  nx2Subset <- floor(nx2*propSubset)
  for (i in 1:nRep) {
    x1subset <- x1[sample(nx1, nx1Subset, replace = FALSE),]
    x2subset <- x2[sample(nx2, nx2Subset, replace = FALSE),]
    x1_mean <- colSums(x1subset)/nx1Subset
    x2_mean <- colSums(x2subset)/nx2Subset
    x1_var <- diag(var(x1subset))*(nx1Subset-1)/nx1Subset
    x2_var <- diag(var(x2subset))*(nx2Subset-1)/nx2Subset
    x_stat_var <- (x1_var*nx1Subset + x2_var*nx2Subset)/(nx1Subset*nx2Subset)
    diff_x_stat <- max((x1_mean/sqrt(x_stat_var) - x2_mean/sqrt(x_stat_var))^2)
    p_x <- 1-exp(-1/sqrt(pi)*exp(-(diff_x_stat-(2*log(p)-log(log(p))))/2))
    if (p_x < 0.05) {
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
  x1Var <- diag(var(x1))*(n1-1)/n1
  x2Var <- diag(var(x2))*(n2-1)/n2
  xStatVar <- (x1Var*n1 + x2Var*n2)/(n1*n2)
  
  SqStandizeDiff <- ((colSums(x1-x2))/n)^2/((diag(var(x1-x2)))*(n-1)/n^2)
  StatX <- max(SqStandizeDiff)
  pvalue <- 1-exp(-1/sqrt(pi)*exp(-(StatX-(2*log(p)-log(log(p))))/2))
  return(pvalue)
}
#------------------------------------------------------------------------------------

require(ggplot2)
require(reshape)

set.seed(150)

# Load the data
setwd("C:/Users/yuanpei/Desktop/two sample test/20160916realdata/dataset/combo")
w <- as.matrix(read.table("bmi_counts.dat"))

# Filter the data
n<-dim(w)[1]
p<-dim(w)[2]

# Add 0.5 on zeros
w[w==0]<-w[w==0]+0.5

# Separate the sample into two groups
BMI_and_Counts <- read.table("C:/Users/yuanpei/Desktop/two sample test/20160916realdata/dataset/combo/BMI_and_Counts.txt", header=T, quote="\"")
bmi<-BMI_and_Counts[,2]

index_1<-which(bmi>=25)
index_2<-which(bmi<25)

wx <- w[index_1,]
wy <- w[index_2,]

nx <- dim(wx)[1]
ny <- dim(wy)[1]

x <- wx/(rowSums(wx)%*%matrix(1,1,p))
y <- wy/(rowSums(wy)%*%matrix(1,1,p))

log_x<-log(x)
log_y<-log(y)

clog_TX <- log_x-1/p*rowSums(log_x)%*%matrix(1,1,p)
clog_TY <- log_y-1/p*rowSums(log_y)%*%matrix(1,1,p)

# Mean
x_mean <- colSums(x)/nx
y_mean <- colSums(y)/ny

log_x_mean <- colSums(log_x)/nx
log_y_mean <- colSums(log_y)/ny

clr_x_mean <- colSums(clog_TX)/nx
clr_y_mean <- colSums(clog_TY)/ny

# Variance
x_var <- diag(var(x))*(nx-1)/nx
y_var <- diag(var(y))*(ny-1)/ny
x_stat_var <- (x_var*nx + y_var*ny)/(nx*ny)

log_x_var <- diag(var(log_x))*(nx-1)/nx
log_y_var <- diag(var(log_y))*(ny-1)/ny
log_x_stat_var <- (log_x_var*nx + log_y_var*ny)/(nx*ny)

clr_x_var <- diag(var(clog_TX))*(nx-1)/nx
clr_y_var <- diag(var(clog_TY))*(ny-1)/ny
clr_x_stat_var <- (clr_x_var*nx + clr_y_var*ny)/(nx*ny)

# p-values
x_stat <- max((x_mean/sqrt(x_stat_var) - y_mean/sqrt(x_stat_var))^2)
p_x <- 1-exp(-1/sqrt(pi)*exp(-(x_stat-(2*log(p)-log(log(p))))/2))
log_x_stat <- max((log_x_mean/sqrt(log_x_stat_var) - log_y_mean/sqrt(log_x_stat_var))^2)
p_logx <- 1-exp(-1/sqrt(pi)*exp(-(log_x_stat-(2*log(p)-log(log(p))))/2))
clr_x_stat <- max((clr_x_mean/sqrt(clr_x_stat_var) - clr_y_mean/sqrt(clr_x_stat_var))^2)
p_clrx <- 1-exp(-1/sqrt(pi)*exp(-(clr_x_stat-(2*log(p)-log(log(p))))/2))

# empirical power
powerX = c()
powerLogX = c()
powerClrX = c()
for (i in 2:10){
  powerX <- c(powerX,powerFun(x,y,i*0.1,5000))
  powerLogX <- c(powerLogX,powerFun(log_x,log_y,i*0.1,5000))
  powerClrX <- c(powerClrX,powerFun(clog_TX,clog_TY,i*0.1,5000))
}

# plot the curve
grid = c(rep(c(2:10)/10, 3))
power = c(powerClrX,powerLogX,powerX)
Method = c(rep(c("Proposed"), 9),rep(c("Log"), 9),rep(c("Raw"), 9))
EmpiricalPower<-data.frame(Grid = grid, Power = power, Method = Method)
setEPS(width=8, height=6)
postscript(file="C:/Users/yuanpei/Desktop/two sample test/20160916realdata/comboPower.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
powerPlot <- qplot(Grid, Power, data = EmpiricalPower, shape = Method, linetype = Method,geom = "path", xlab = "Proportion of sample size", ylab = "Power")
powerPlot + scale_colour_grey() + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.title = element_text(face="bold", 
                                                                                                                                                     size="12"), legend.position="right")+geom_point(size=3)
dev.off()
