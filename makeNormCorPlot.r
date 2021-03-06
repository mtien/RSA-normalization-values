
rm(list = ls())
graphics.off()

options(scipen=-1)

library(ggplot2)
require(grid)
library(plyr)

normalize <- function(vect) {
	(vect-min(vect))/(max(vect)-min(vect))
}


multiplot <- function(..., plotlist=NULL, cols) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  print(numPlots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}

plotter <- function(dat, y.label = 'H-Scale', x.label = 'Average RSA', 
                     op.title = '', pval = 0, cor = 0) {
  
  max_figure_width_si <- 10
  the_pointsize <- 17.5
  

  theme_set(theme_bw(base_size=the_pointsize))
  old_theme <- theme_update(panel.border=theme_blank(),
                            axis.line=theme_segment(),
                            panel.grid.minor=theme_blank(),
                            panel.grid.major=theme_blank(),
                            panel.background=theme_blank(),
                            panel.border=theme_blank(),
                            axis.line=theme_segment())
  
  g <- ggplot(dat, aes(x=x, y=y))
  g <- g + geom_point(col='black', size=2.5)
  g <- g + opts(strip.background=theme_blank())
  g <- g + opts(panel.border=theme_blank(), axis.line=theme_segment())
  g <- g + ylab(y.label)
  g <- g + xlab(x.label)
  g <- g + scale_x_continuous(breaks=c(0, .5, 1), limits=c(0, 1.15))
  g <- g + scale_y_continuous(breaks=c(0,.25, .5,.75, 1), limits=c(0,1.15))
  g <- g + opts(axis.title.x = theme_text(size=the_pointsize, vjust=-.5))
  g <- g + opts(axis.title.y = theme_text(angle=90, size=the_pointsize, vjust=.3))
  g <- g + opts(title=op.title)
  g <- g + opts(plot.margin=unit(c(0, 1, 1, 1), "lines"))
  g <- g + opts(axis.ticks.margin = unit(0.15, "cm"))
  g <- g + opts(legend.position = "none")
  g <- g + geom_text( x=.5, y=1.15, label= paste('r = ',round(cor, digits=3),
                                                              sep=''))
  return(g)
}

setwd("~/GitHub/RSA-normalization-values/Correlation/")
wolfenden<- read.delim("Wolfden.txt")
wolf<- normalize(wolfenden$Hydro)

Rose<-read.delim("rose.txt")
rose<- normalize(Rose$Hydro)

Kyte<-read.delim("Kite_Doolittle.txt")
kyte<- normalize(Kyte$Hydro)

Scales<- read.table("Hydrophobicity_Scales_updated.txt", row.names=1, header=TRUE)
meanTheo<-normalize(Scales$meanTheo)
medianTheo<-normalize(Scales$medianTheo)
SqrtMeanTheo<-normalize(Scales$meanTheoSqrt)
Frac<- normalize(Scales$ScaleFrac)

meanEmp<-normalize(Scales$meanEmp)
medianEmp<-normalize(Scales$medianEmp)
SqrtMeanEmp<-normalize(Scales$meanEmpSqrt)

Wolf.TheoMean<-data.frame( x= meanTheo, y= wolf) 
Wolf.Rose<-data.frame( x= rose, y= wolf)
Wolf.TheoMedian<-data.frame( x= medianTheo, y= wolf)
Wolf.TheoSqrt<-data.frame( x= SqrtMeanTheo, y= wolf)
Wolf.Fraction<-data.frame( x= Frac, y= wolf)

Wolf.EmpMean<-data.frame( x= meanEmp, y= wolf) 
Wolf.EmpMedian<-data.frame( x= medianEmp, y= wolf)
Wolf.EmpSqrt<-data.frame( x= SqrtMeanEmp, y= wolf)

Kyte.TheoMean<-data.frame( x= meanTheo, y=kyte)
Kyte.Rose<-data.frame( x= rose, y=kyte)
Kyte.TheoMedian<-data.frame( x= medianTheo, y=kyte)
Kyte.TheoSqrt<-data.frame( x= SqrtMeanTheo, y=kyte)
Kyte.Fraction<-data.frame( x= Frac, y=kyte)

Kyte.EmpMean<-data.frame( x= meanEmp, y=kyte)
Kyte.EmpMedian<-data.frame( x= medianEmp, y=kyte)
Kyte.EmpSqrt<-data.frame( x= SqrtMeanEmp, y=kyte)

cor.wolf.theoMean<-cor.test(wolf, meanTheo)
cor.wolf.theoMedian<-cor.test(wolf, medianTheo)
cor.wolf.theoSqrt<-cor.test(wolf, SqrtMeanTheo)
cor.wolf.rose<-cor.test(wolf, rose)
cor.wolf.fraction<-cor.test(wolf, Frac)

cor.wolf.empMean<-cor.test(wolf, meanEmp)
cor.wolf.empMedian<-cor.test(wolf, medianEmp)
cor.wolf.empSqrt<-cor.test(wolf, SqrtMeanEmp)

cor.kyte.theoMean<-cor.test(kyte, meanTheo)
cor.kyte.theoMedian<-cor.test(kyte, medianTheo)
cor.kyte.theoSqrt<-cor.test(kyte, SqrtMeanTheo)
cor.kyte.rose<-cor.test(kyte, rose)
cor.kyte.fraction<-cor.test(kyte, Frac)

cor.kyte.empMean<-cor.test(kyte, meanEmp, method='sp')
cor.kyte.empMedian<-cor.test(kyte, medianEmp, method='sp')
cor.kyte.empSqrt<-cor.test(kyte, SqrtMeanEmp, method='sp')

p.Wolf.TheoMean<-plotter(Wolf.TheoMean, 
					y.label='Wolfenden H-scale', 
					x.label='Mean RSA', 
					cor= cor.wolf.theoMean$estimate,
					pval= cor.wolf.theoMean$p.value)
##p.Wolf.TheoMedian<-plotter(Wolf.TheoMedian, 
##					y.label='Wolfenden H-scale', 
##					x.label='Median RSA', 
##					cor= cor.wolf.theoMedian$estimate,
##					pval= cor.wolf.theoMedian$p.value) 
##p.Wolf.TheoSqrt<-plotter(Wolf.TheoSqrt, 
##					y.label='Wolfenden H-scale', 
##					x.label='Mean Squareroot RSA', 
##					cor= cor.wolf.theoSqrt$estimate,
##					pval= cor.wolf.theoSqrt$p.value)
p.Wolf.Rose<-plotter(Wolf.Rose, 
					y.label='Wolfenden H-scale', 
					x.label='Rose H-scale', 
					cor= cor.wolf.rose$estimate,
					pval= cor.wolf.rose$p.value)
p.Wolf.Fraction<-plotter(Wolf.Fraction, 
					y.label='Wolfenden H-scale', 
					x.label='Fraction of 100% Buried Residues', 
					cor= cor.wolf.fraction$estimate,
					pval= cor.wolf.fraction$p.value)
##p.Wolf.EmpMean<-plotter(Wolf.EmpMean, 
##					y.label='Wolfenden H-scale', 
##					x.label='Mean RSA', 
##					cor= cor.wolf.empMean$estimate,
##					pval= cor.wolf.empMean$p.value)
##p.Wolf.EmpMedian<-plotter(Wolf.EmpMedian, 
##					y.label='Wolfenden H-scale', 
##					x.label='Median RSA', 
##					cor= cor.wolf.theoMedian$estimate,
##					pval= cor.wolf.theoMedian$p.value) 
##p.Wolf.TheoSqrt<-plotter(Wolf.TheoSqrt, 
##					y.label='Wolfenden H-scale', 
##					x.label='Mean Squareroot RSA', 
##					cor= cor.wolf.theoSqrt$estimate,
##					pval= cor.wolf.theoSqrt$p.value)
##
##
p.Kyte.TheoMean<-plotter(Kyte.TheoMean, 
					y.label='Kyte-Doolittle H-scale', 
					x.label='Mean RSA', 
					cor= cor.kyte.theoMean$estimate,
					pval= cor.kyte.theoMean$p.value)
##p.Kyte.TheoMedian<-plotter(Kyte.TheoMedian, 
##					y.label='Kyte H-scale', 
##					x.label='Median RSA', 
##					cor= cor.kyte.theoMedian$estimate,
##					pval= cor.kyte.theoMedian$p.value) 
##p.Kyte.TheoSqrt<-plotter(Kyte.TheoSqrt, 
##					y.label='Kyte H-scale', 
##					x.label='Mean Squareroot RSA', 
##					cor= cor.kyte.theoSqrt$estimate,
##					pval= cor.kyte.theoSqrt$p.value)
p.Kyte.Rose<-plotter(Kyte.Rose, 
					y.label='Kyte-Doolittle H-scale', 
					x.label='Rose H-scale', 
					cor= cor.kyte.rose$estimate,
					pval= cor.kyte.rose$p.value)
p.Kyte.Fraction<-plotter(Kyte.Fraction, 
					y.label='Kyte-Doolittle H-scale', 
					x.label='Fraction of 100% Buried Residues', 
					cor= cor.kyte.fraction$estimate,
					pval= cor.kyte.fraction$p.value)
##p.Kyte.EmpMean<-plotter(Kyte.EmpMean, 
##					y.label='Kyte H-scale', 
##					x.label='Mean RSA', 
##					cor= cor.kyte.empMean$estimate,
##					pval= cor.kyte.empMean$p.value)
##p.Kyte.EmpMedian<-plotter(Kyte.EmpMedian, 
##					y.label='Kyte H-scale', 
##					x.label='Median RSA', 
##					cor= cor.kyte.theoMedian$estimate,
##					pval= cor.kyte.theoMedian$p.value) 
##p.Kyte.EmpSqrt<-plotter(Kyte.TheoSqrt, 
##					y.label='Kyte H-scale', 
##					x.label='Mean Squareroot RSA', 
##					cor= cor.kyte.theoSqrt$estimate,
##					pval= cor.kyte.theoSqrt$p.value)

setwd("~/GitHub/RSA-normalization-values/Figures/")
pdf('NormalizedCorrelations.pdf', width=15, height=8 ,useDingbats=F)

multiplot(p.Wolf.Rose, p.Wolf.TheoMean,  p.Wolf.Fraction,
          p.Kyte.Rose, p.Kyte.TheoMean,  p.Kyte.Fraction,
          cols=3)

dev.off()
