library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)

#read all the files
flist<-c("ctpgradient.txt")
dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data<-data.frame(dt_list)
names(data)<- c("psize", "eps010", "eps011", "eps012", "eps013", "eps014")

df <- melt(data, id.vars="psize")


myplot <-ggplot(df, aes(x=psize, y=value, colour=variable))+
  geom_line(size=2)
myplot+
  geom_hline(yintercept = 0.0, colour="black", linetype="dashed")+
  geom_vline(xintercept = 250, colour="black", linetype="dashed")+
  xlab("Number of (5,3) strategies (pop. size=500)")+
  ylab("Gradient of selection (T+ - T-)")+
  scale_color_manual(name="Reasoning noise",values=c("#D55E00","#E69F00","#009E73","#0072B2","#000000"),
                        labels=c(epsilon~"=0.05", epsilon~"=0.06", epsilon~"=0.07", epsilon~"=0.08", epsilon~"=0.09"))+
  theme_bw()+
  theme(legend.position = c(0.2, 0.8),
        legend.text = element_text(size=25),
        legend.title = element_text(size=25),
        axis.title.x=element_text(size=25),
        axis.title.y=element_text(size=25),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_text(size=25),
        axis.text.y=element_text(size=25))
  
