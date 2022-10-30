library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggtext)

#read all the files

flist<-c("ctpbeliefsS10.txt","ctpbeliefsS01.txt","ctpbeliefsS001.txt")

dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data<-data.frame(do.call(rbind,dt_list))
names(data)<-c("Scale","1","2","3","4","5")
data[1,1]<-"1"
data[2,1]<-"10^-1"
data[3,1]<-"10^-2"


df <- melt(data, id.vars="Scale")
df[,3] <- round (df[,3],2)

plot<-ggplot(data=df, aes(x=variable, y=value, fill=factor(Scale))) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  labs(x="belief(*t*)", y = "Frequency", fill="scale")+
  geom_text(aes(label=value), vjust=1.6, color="white",
            position = position_dodge(0.9), size=7)+
  scale_fill_brewer(labels=c("*\u03BB=1*",  "*\u03BB=10^-1*", "*\u03BB=10^-2*"), palette="Blues", direction=1)+
  theme_minimal()+
  guides(fill = guide_legend(byrow = TRUE))+
  theme(legend.position = c(0.2, 0.8),
        legend.text = element_markdown(size=27),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.5,"cm"),
        axis.title.x = element_markdown(size=35),
        axis.title.y=element_text(size=35),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_text(size=27),
        axis.text.y=element_text(size=27))
plot

