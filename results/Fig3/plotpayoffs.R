library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)

#read all the files
flist<-c("ctppayoffsLall.txt")
dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data<-data.frame(dt_list)
names(data)<- c("eps", "(t=1,k=0)", "(t=5,k=3)")
df <- melt(data, id.vars="eps")


myplot <-ggplot(df, aes(x=eps, y=value, colour=variable, fill=variable))+
  geom_line(aes(group=variable))+
  geom_point(size=4, shape=21)
myplot+coord_fixed()+
  geom_vline(xintercept = 0.13, colour="black", linetype="dashed")+
  annotate("text", x=0.131, y=0.73, label= "("~epsilon~"=0.13)", colour="black", size=5, hjust="left")+
  xlab(epsilon~"(reasoning noise)")+
  ylab("Fitness")+
  scale_fill_manual(name="Strategy",values=c("#D55E00","#000000"),
                     labels=c("(t=1,k=0)", "(t=5,k=3)"))+
  scale_color_manual(name="Strategy",values=c("#D55E00","#000000"),
                    labels=c("(t=1,k=0)", "(t=5,k=3)"))+
  theme_bw()+
  theme(legend.position = c(0.23, 0.8),
        legend.text = element_text(size=25),
        legend.title = element_text(size=25),
        axis.title.x=element_text(size=25),
        axis.title.y=element_text(size=25),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_text(size=25),
        axis.text.y=element_text(size=25))
  
