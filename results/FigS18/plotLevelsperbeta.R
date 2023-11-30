library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)

#read all the files
flist<-list.files(pattern="ctplevelsperbeta.txt")
dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data<-data.frame(dt_list)
names(data)<- c("beta", "level","0", "1", "2", "3", "4", "5", "6", "total")
datasub<-select(data, beta, level,total)


myplot <-ggplot(datasub, aes(x=beta, y=total, group=as.factor(level)))+
  geom_line(aes(color=as.factor(level)))+
  geom_point(size=5, aes(shape=as.factor(level),color=as.factor(level),fill=as.factor(level)))+
  scale_color_manual(name="Level (k)", 
                     values=c("#D55E00","#E69F00", "#009E73","#0072B2","#000000","#56B4E9","#CC79A7"),
                     labels=c("0", "1", "2", "3", "4","5","6"))+
  scale_fill_manual(name="Level (k)", 
                    values=c("#D55E00","#E69F00", "#009E73","#0072B2","#000000","#56B4E9","#CC79A7"),
                    labels=c("0", "1", "2", "3", "4","5","6"))+
  scale_shape_manual(name="Level (k)", 
                     values=c(21, 22, 23,24,25, 0, 1),
                     labels=c("0", "1", "2", "3", "4","5","6"))+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  geom_vline(xintercept = 0.04, colour="black", linetype="dashed")+
  annotation_logticks() +
  xlab(beta~"(selection strength)")+
  coord_cartesian(xlim=c(0.0001,10.0))+
  ylab("Frequency")+
  theme_bw()+
  theme(legend.position = c(0.2, 0.8),
        legend.text = element_text(size=27),
        legend.title = element_text(size=27),
        axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_text(size=27),
        axis.text.y=element_text(size=27))
myplot

