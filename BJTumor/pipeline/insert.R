#!/usr/bin/Rscript
#xinghe <xingh3223@berryoncology.com>

library(ggplot2)
library(scales)
args = commandArgs(trailingOnly=TRUE)
data <- read.table(args[1], header=FALSE)
png(args[2])
data <- as.data.frame(data)
#plot(t(data),type='h', xlim=c(0,400), ylim=c(0,max(data[1,])), xlab='Postion of Reads', ylab='Error rate (%)', col='#34659F', xaxs='i', yaxs='i',las=1,main="Distribution of base error rate")

#plot(data,type='l', xlim=c(0,max(data[,1])), ylim=c(0,max(data[,2])), xlab='Insert size (bp)', ylab='Percent', col='#34659F', xaxs='i', yaxs='i',las=1)
p <- ggplot(data,aes(x=V1,y=V2))+geom_line(color='#34659F')+scale_y_continuous(labels=percent_format())+labs(x='Insert Size (bp)',y='Percent')
p <-p + theme_light()+theme(panel.grid.major=element_blank())
p
dev.off()
