grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}
#-----------------------------------------------------------------------------------------------------------------

require(ggplot2)
require(reshape)
require(xlsx)
require(gridExtra)
require(grid)
require(cowplot)
require(Cairo)
require(latex2exp)

setHook(packageEvent("grDevices", "onLoad"),
        function(...) grDevices::X11.options(type='cairo'))
options(device='x11')

#----------------------------------Power Curve Plot---------------------------------------------------------------
# Load the table
# Load Combo power curve
xt1 <- as.matrix(read.xlsx("table/power-bootstrap-backtest.xlsx", sheetName="Combo-power"))
powerXC <- as.double(xt1[3,-1])
powerLogXC <- as.double(xt1[2,-1])
powerClrXC <- as.double(xt1[1,-1])

grid <- c(rep(c(4:20)/20, 3)) # set the grid 0.2 : 0.05 : 1
powerC <- c(powerClrXC,powerLogXC,powerXC)
MethodC <- c(rep(c("Proposed"), 17),rep(c("Log"), 17),rep(c("Raw"), 17))
EmpiricalPowerC <- data.frame(Grid = grid, Power = powerC, Method = MethodC)

# Load Please power curve
xt2 <- as.matrix(read.xlsx("table/power-bootstrap-backtest.xlsx", sheetName="Please-power-p43-1vs4"))
powerXP <- as.double(xt2[3,-1])
powerLogXP <- as.double(xt2[2,-1])
powerClrXP <- as.double(xt2[1,-1])

grid = c(rep(c(8:20)/20, 3))
powerP = c(powerClrXP,powerLogXP,powerXP)
MethodP = c(rep(c("Proposed"), 13),rep(c("Log"), 13),rep(c("Raw"), 13))
EmpiricalPowerP <- data.frame(Grid = grid, Power = powerP, Method = MethodP)

# Load back-testing histogram
bkTestC <- data.frame(read.xlsx("table/power-bootstrap-backtest.xlsx", sheetName = "Combo-backtesting", header = FALSE))
bkTestP <- data.frame(read.xlsx("table/power-bootstrap-backtest.xlsx", sheetName = "Please-backtesting-p43-1vs4", header = FALSE))


# Plot two power curves in one figure
#setEPS(width = 10, height = 7)
#postscript(file="Plot/comboPleaseFig.eps", family="Times")
#cairo_pdf(file="Plot/comboPleaseFig.pdf", family="Times")


setEPS(width = 4.25, height = 3.5)
postscript(file="Plot/comboPleaseFig1.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
powerPlot1 <- ggplot(data = EmpiricalPowerC, aes(x = Grid, y = Power, shape = Method, linetype = Method)) + 
  geom_path() + geom_point() +  xlab("Subsampling proportion") + ylab("Power") + 
  scale_colour_grey() + theme_bw() + theme(panel.border = element_rect(size = 1, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
#   ggtitle('(a)') + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),axis.title = element_text(size = 15), plot.title = element_text(face = "plain", hjust = 0)) 
#  scale_x_continuous(breaks = c(0.2,0.4,0.6,0.8,1.0), labels = c(TeX('$0\\cdot 20$'),TeX('$0\\cdot 40$'),TeX('$0\\cdot 60$'),TeX('$0\\cdot 80$'),TeX('$1\\cdot 00$')))  + 
#  scale_y_continuous(limits = c(0.1, 1), breaks = c(0.25,0.5,0.75,1.00), labels = c(TeX('$0\\cdot 25$'),TeX('$0\\cdot 50$'),TeX('$0\\cdot 75$'),TeX('$1\\cdot 00$'))) 
powerPlot1
dev.off()

setEPS(width = 4.25, height = 3.5)
postscript(file="Plot/comboPleaseFig2.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
powerPlot2 <- ggplot(data = EmpiricalPowerP, aes(x = Grid, y = Power, shape = Method, linetype = Method)) + 
  geom_path() + geom_point() +  xlab("Subsampling proportion") + ylab("Power") + 
  scale_colour_grey() + theme_bw() + theme(panel.border = element_rect(size = 1, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
#   ggtitle('(b)') + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),axis.title = element_text(size = 15), plot.title = element_text(face = "plain", hjust = 0)) 
#  scale_x_continuous(breaks = c(0.4,0.6,0.8,1.0), labels = c(TeX('$0\\cdot 40$'),TeX('$0\\cdot 60$'),TeX('$0\\cdot 80$'),TeX('$1\\cdot 00$'))) +
#  scale_y_continuous(limits = c(0.1, 1), breaks = c(0.25,0.5,0.75,1.00), labels = c(TeX('$0\\cdot 25$'),TeX('$0\\cdot 50$'),TeX('$0\\cdot 75$'),TeX('$1\\cdot 00$'))) 
powerPlot2
dev.off()


bkTestCPlot <- ggplot(bkTestC, aes(X1)) + geom_histogram(fill = "white", color = "black", bins = 10) + ylim(0, 200) + xlab(expression(italic(p)*"-"*value)) + 
  ylab("Frequency") +  
  # ggtitle('(c)') + 
  theme(panel.border = element_rect(size = 1, color = "black"), panel.grid.major = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),axis.title = element_text(size = 15), plot.title = element_text(face = "plain", hjust = 0)) 
#  scale_x_continuous(breaks = c(0.0,0.25,0.5,0.75,1.00), labels = c(TeX('$0\\cdot 00$'),TeX('$0\\cdot 25$'),TeX('$0\\cdot 50$'),TeX('$0\\cdot 75$'),TeX('$1\\cdot 00$'))) 

bkTestPPlot <- ggplot(bkTestP, aes(X1)) + geom_histogram(fill = "white", color = "black", bins = 10) + ylim(0, 200) + xlab(expression(italic(p)*"-"*value)) + 
  ylab("Frequency") +  
  #ggtitle('(d)') + 
  theme(panel.border = element_rect(size = 1, color = "black"), panel.grid.major = element_blank(), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title = element_text(size = 15), plot.title = element_text(face = "plain", hjust = 0))  
#  scale_x_continuous(breaks = c(0.0,0.25,0.5,0.75,1.00), labels = c(TeX('$0\\cdot 00$'),TeX('$0\\cdot 25$'),TeX('$0\\cdot 50$'),TeX('$0\\cdot 75$'),TeX('$1\\cdot 00$'))) 

grid.arrange(powerPlot1, powerPlot2, bkTestCPlot, bkTestPPlot, ncol = 2)
#grid_arrange_shared_legend(powerPlot1, powerPlot2, ncol=2)
dev.off()

#----------------------------------Bootstrapping and Backtesing Histogram------------------------------------------
# Load the table
# Load Combo boostrapping and backtesting histogram
# bsC <- as.matrix(read.xlsx("table/power-bootstrap-backtest.xlsx", sheetName = "Combo-bootsraping", header = FALSE))
bkTestC <- as.matrix(read.xlsx("table/power-bootstrap-backtest.xlsx", sheetName = "Combo-backtesting", header = FALSE))
# bsP <- as.matrix(read.xlsx("table/power-bootstrap-backtest.xlsx", sheetName = "Please-bootstraping-p43-1vs4", header = FALSE))
bkTestP <- as.matrix(read.xlsx("table/power-bootstrap-backtest.xlsx", sheetName = "Please-backtesting-p43-1vs4", header = FALSE))


#library(Cairo)
#cairo_ps(filename = 'Plot/comboPleaseHistogram.ps', width = 9, height = 4, onefile = T, family = 'Times')
#cairo_pdf(filename = 'Plot/comboPleaseHistogram.pdf', width = 9, height = 4, onefile = T, family = 'Times')
#par(mfrow = c(1,2))
#hist(bsC, xlim=c(0,1.0), ylim=c(0,1000), xlab = "p-value", ylab = "Frequency", main = NULL, col = rgb(0.1,0.1,0.1,0.5), cex.lab = 1.25)

setEPS(width = 4.25, height = 3.5)
postscript(file="Plot/comboPleaseFig3.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
hist(bkTestC, main = NULL, xlim = c(0.0,1.0), ylim=c(0,250), xlab = expression(italic(p)*"-"*value), ylab = "Frequency", col = rgb(0.8,0.8,0.8,0.5), cex.lab = 1.25)
#legend(pt.cex = 2, pch = 15)
#hist(bsP, xlim=c(0,1.0), ylim=c(0,1000), xlab = "p-value", ylab = "Frequency", main = NULL, col = rgb(0.1,0.1,0.1,0.5), cex.lab = 1.25)
dev.off()

setEPS(width = 4.25, height = 3.5)
postscript(file="Plot/comboPleaseFig4.eps", family="Times")
par(mai=c(0.6, 0.6, 0.1, 0), mgp=c(1.8, 0.6, 0))
hist(bkTestP, main = NULL, xlim = c(0.0,1.0), ylim=c(0,200), xlab = expression(italic(p)*"-"*value), ylab = "Frequency", col = rgb(0.8,0.8,0.8,0.5), cex.lab = 1.25)
#legend(pt.cex = 2, pch = 15)
#dev.off()
dev.off()

plottmp <- ggplot(df, aes(X1)) + geom_histogram(fill = "white", color = "black", bins = 10) + ylim(0, 200) + xlab("p-value") + ylab("Frequency")
