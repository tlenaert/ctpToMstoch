library(magrittr)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggtext)


#read all the files
flist<-c("ctpdecperKLall125.txt","ctpdecperKLall250.txt","ctpdecperKLall500.txt")

dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting

data<-data.frame(do.call(rbind,dt_list))
names(data)<-c("psize","level","1","2","3","4","5","total")
data[1:5,1]<-"125"
data[6:10,1]<-"250"
data[11:15,1]<-"500"
datasub <-select(data, psize, level,"1","2","3","4","5")

df <- melt(datasub, id.vars=c("level","psize"))
df[,4] <- round (df[,4],2)

plot<-ggplot(data=df, aes(x=fct_inorder(factor(psize)), y=value, fill=factor(variable))) +
  geom_bar(stat="identity", color="black")+
  labs(title="*k*-levels", x="Population size *Z*", y = "Frequency", fill="actions(*T*)")+
  scale_fill_brewer(palette="Oranges")+
  facet_grid(~level)+
  theme_minimal()+
  guides(fill = guide_legend(byrow = TRUE))+
  theme(legend.position = c(0.25, 0.8),
        legend.text = element_markdown(size=27),
        legend.title = element_markdown(size=30),
        legend.spacing.y = unit(0.1,"cm"),
        plot.title = element_markdown(size=35, hjust = 0.5),
        axis.title.x=element_markdown(size=35),
        axis.title.y=element_text(size=35),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_markdown(size=20),
        axis.text.y=element_markdown(size=27),
        strip.text.x=element_text(size=27))
plot

 