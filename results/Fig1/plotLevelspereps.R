library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)

#read all the files
flist<-list.files(pattern="ctplevelspereps.txt")
dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data<-data.frame(dt_list)
names(data)<- c("eps", "level","0", "1", "2", "3", "4", "total")
datasub<-select(data, eps, level,total)


myplot <-ggplot(datasub, aes(x=eps, y=total, group=as.factor(level)))+
  geom_line(aes(color=as.factor(level)))+
  geom_point(size=5, aes(shape=as.factor(level),color=as.factor(level),fill=as.factor(level)))+
  scale_color_manual(name="Level (k)", 
                     values=c("#D55E00","#E69F00","#009E73","#0072B2","#000000"),
                     labels=c("0", "1", "2", "3", "4"))+
  scale_fill_manual(name="Level (k)", 
                    values=c("#D55E00","#E69F00","#009E73","#0072B2","#000000"),
                    labels=c("0", "1", "2", "3", "4"))+
  scale_shape_manual(name="Level (k)", 
                    values=c(21, 22, 23,24,25),
                    labels=c("0", "1", "2", "3", "4"))+
  expand_limits(x=0.5)+
  geom_vline(xintercept = 0.18, colour="black", linetype="dashed")+
  xlab(epsilon~"(reasoning error probability)")+
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


