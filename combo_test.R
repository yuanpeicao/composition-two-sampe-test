require(ggplot2)
require(reshape)
# Load the data
#setwd("C:/Users/yuanpei/Desktop/two sample test/20160916realdata/dataset")
w <- as.matrix(read.table("~/comptest_code/comptest_realdata/dataset/combo/bmi_counts.dat"))

# Filter the data
n<-dim(w)[1]
#binary_w<-w
#binary_w[binary_w>0]<-1
#data_filter<-5/n
#index_w<-which(colMeans(binary_w)>data_filter)
#w<-w[,index_w]
#w<-w[,colMeans(w)>data_filter]

# # Add 0.5 on zeros
# w[w==0]<-w[w==0]+0.5

# Add 0.5 on zeros
w[w==0]<-w[w==0]+0.1

# Delete the 58th bacteria (Allisonella)
#w <- w[,-58]
#w <- w[,-36]
p<-dim(w)[2]

# Separate the sample into two groups
quan<-0.25
BMI_and_Counts <- read.table("~/comptest_code/comptest_realdata/dataset/combo/BMI_and_Counts.txt", header=T, quote="\"")
bmi<-BMI_and_Counts[,2]

#index_1<-which(bmi>quantile(bmi,1-quan))
#index_2<-which(bmi<=quantile(bmi,quan))

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

## Generate barplot
# difference on X
high_BMI_x <- x_mean/sqrt(x_stat_var)
low_BMI_y <- y_mean/sqrt(x_stat_var)
x_stat <- max((high_BMI_x - low_BMI_y)^2)
x_ind <- which((high_BMI_x - low_BMI_y)^2 == x_stat)
p_x <- 1-exp(-1/sqrt(pi)*exp(-(x_stat-(2*log(p)-log(log(p))))/2))

#species_index<-as.matrix(read.table("bmi_species_filter_5.dat"))
species_index <- as.matrix(read.table("~/comptest_code/comptest_realdata/dataset/combo/bmi_species_whole.dat"))
Genera <- paste(species_index,sep="")
#Genera <- Genera[-58]
#Genera <- Genera[-36]

BMI_data_x <- data.frame(Genera=Genera,"Obese group"=high_BMI_x,"Lean group"=low_BMI_y)
melted_BMI_data_x <- melt(BMI_data_x, id.vars=c("Genera"))
ggplot(melted_BMI_data_x, aes(Genera, value, fill=variable)) + 
  geom_bar(stat="identity", position="dodge")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(y = "standardized X")+ theme(legend.title=element_blank())

# difference on log(X)
high_BMI_logx <- log_x_mean/sqrt(log_x_stat_var)
low_BMI_logy <- log_y_mean/sqrt(log_x_stat_var)
log_x_stat <- max((high_BMI_logx - low_BMI_logy)^2)
log_x_ind <- which((high_BMI_logx - low_BMI_logy)^2 == log_x_stat)
p_log_x <- 1-exp(-1/sqrt(pi)*exp(-(log_x_stat-(2*log(p)-log(log(p))))/2))

BMI_data_logx <- data.frame(Genera=Genera,"Obese group"=high_BMI_logx,"Lean group"=low_BMI_logy)
melted_BMI_data_logx <- melt(BMI_data_logx, id.vars=c("Genera"))
ggplot(melted_BMI_data_logx, aes(Genera, value, fill=variable)) + 
  geom_bar(stat="identity", position="dodge")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(y = "standardized logX")+ theme(legend.title=element_blank())


# difference on clr(x)
high_BMI_clr_x <- clr_x_mean/sqrt(clr_x_stat_var)
low_BMI_clr_y <- clr_y_mean/sqrt(clr_x_stat_var)
diff_clr <- (clr_x_mean - clr_y_mean)^2
clr_x_stat <- max((high_BMI_clr_x - low_BMI_clr_y)^2)
p_clr_x <- 1-exp(-1/sqrt(pi)*exp(-(clr_x_stat-(2*log(p)-log(log(p))))/2))

clr_x_ind <- which((high_BMI_clr_x - low_BMI_clr_y)^2 == clr_x_stat)
clr_bacteria <- Genera[clr_x_ind]
#BMI_data_clrx<-data.frame(Genera=Genera,"Obese group"=high_BMI_clr_x,"Lean group"=low_BMI_clr_y)
BMI_data_clrx<-data.frame(Genera=Genera,clrX=diff_clr)
# BMI_data_clrx<-data.frame(Genera=Genera,clrX=(high_BMI_clr_x-low_BMI_clr_y)^2)
melted_BMI_data_clrx <- melt(BMI_data_clrx, id.vars=c("Genera"))
setEPS()
postscript("barplot_clr_p86.eps", horizontal = TRUE, onefile = TRUE, height = 5, width = 11)
ggplot2(melted_BMI_data_clrx + geom_bar(stat="identity", position="dodge")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(y = "standardized clrX")+ theme(legend.title=element_blank())

ggplot2.barplot(data=BMI_data_clrx, xName="time", yName='total_bill', orientation="yAxisReversed")

dev.off()

#ggplot(melted_BMI_data_clrx, aes(Genera, value, fill=variable)) + 
#  geom_bar(stat="identity", position="dodge")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(y = "standardized clrX")+ theme(legend.title=element_blank())

