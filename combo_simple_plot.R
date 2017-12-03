require(ggplot2)
require(reshape)
# Load the data
setwd("C:/Users/yuanpei/Desktop/two sample test/20160916realdata/dataset")
w <- as.matrix(read.table("bmi_counts.dat"))

# Filter the data
n<-dim(w)[1]

# Add 0.5 on zeros
w[w==0]<-w[w==0]+0.5

# Delete the 58th bacteria (Allisonella)
#w <- w[,-58]
#w <- w[,-36]
p<-dim(w)[2]

# Separate the sample into two groups
quan<-0.25
BMI_and_Counts <- read.table("C:/Users/yuanpei/Desktop/two sample test/20160916realdata/dataset/BMI_and_Counts.txt", header=T, quote="\"")
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

#species_index<-as.matrix(read.table("bmi_species_filter_5.dat"))
species_index <- as.matrix(read.table("bmi_species_whole.dat"))
Genera <- paste(species_index,sep="")
genera.info <-as.matrix(read.csv("genera87.csv"))
taxa <- matrix(unlist(strsplit(colnames(genera.info), "\\.")), 6)
phyla <- taxa[2, ]
#Genera <- Genera[-58]
#Genera <- Genera[-36]

# difference on clr(x)
obese_clrx <- clr_x_mean/sqrt(clr_x_stat_var)
lean_clrx <- clr_y_mean/sqrt(clr_x_stat_var)
diff_clrx <- (obese_clrx - lean_clrx)^2
clr_x_stat <- max((obese_clrx - lean_clrx)^2)
p_clrx <- 1-exp(-1/sqrt(pi)*exp(-(clr_x_stat-(2*log(p)-log(log(p))))/2))
clr_x_ind <- which((obese_clrx - lean_clrx)^2 == clr_x_stat)
clr_bacteria <- Genera[clr_x_ind]

BMI_data_clrx <- data.frame(Genera = Genera, clrX = diff_clrx, Phylum = phyla)
BMI_data_clrx <- BMI_data_clrx[order(BMI_data_clrx$clrX,decreasing = TRUE),]
BMI_data_clrx <- BMI_data_clrx[1:10,]
melted_BMI_clrx <- melt(BMI_data_clrx, id.vars=c("Genera"))

setEPS(width=8, height=4)
postscript(file="C:/Users/yuanpei/Desktop/two sample test/20160916realdata/combo.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
p <- ggplot(BMI_data_clrx, aes(x=reorder(Genera, clrX),y=clrX))
p <- p + geom_bar(stat="identity", position = "identity", width=0.5)
p <- p + labs(x="", y="Squared standardized difference on Y")
p <- p + scale_y_continuous(breaks=seq(0,25,by=2))+ theme(panel.grid.major=element_line(colour=NA), text = element_text(size=15), axis.text.x = element_text(size  = 15, angle = 90, hjust = 1, vjust = 1))
p <- p + coord_flip()
p
dev.off()

#ggplot(BMI_data_clrx, aes(Genera, clrX, fill=Phylum)) + 
#  geom_bar(stat="identity", position="dodge")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(y = "Standardized difference")+ theme(legend.title=element_blank()) + coord_flip()
