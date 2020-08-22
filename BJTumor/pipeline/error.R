#!/usr/bin/Rscript
#xinghe <xingh3223@berryoncology.com>
library(ggplot2)
library(scales)
args = commandArgs(trailingOnly=TRUE)
data <- read.table(args[1], header=FALSE, row.names=1)
png(args[2])
data=t(data)
data <-as.data.frame(data)
#plot(data,type='h', xlim=c(0,max(data[,1])), ylim=c(0,max(data[,2])), xlab='Postion of Reads', ylab='Error rate (%)', col='#34659F', xaxs='i', yaxs='i',las=1,main="Distribution of base error rate")
p <-ggplot(data=data,aes(x=pos,y=error))
p <-p +geom_area(fill='#34659F')+geom_line(color='#34659F')+scale_y_continuous(labels=percent_format())+labs(x='Postion of Reads',y='Error Rate')
p <- p+theme_light()+theme(panel.grid.major=element_blank())
p
dev.off()
