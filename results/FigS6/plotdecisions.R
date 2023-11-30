library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggtext)

#read all the files
flist<-list.files(pattern="ctpdecperKLall.txt")

dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data<-data.frame(dt_list)
names(data)<- c("beta", "level","1", "2", "3", "4","5","6","7","total")
datasub<-select(data,"level","1", "2","3","4","5","6","7")
df <- melt(datasub, id.vars="level")
df[,3] <- round (df[,3],2)


plot<-ggplot(data=df, aes(x=level, y=value, fill=variable)) +
  geom_bar(stat="identity", color="black")+
  labs(x="*k*-levels", y = "Frequency", fill="actions(*T*)")+
  scale_fill_brewer(palette="Oranges")+
  scale_x_continuous(breaks=seq(0,6,1))+
  theme_minimal()+
  guides(fill = guide_legend(byrow = TRUE))+
  theme(legend.position = c(0.65, 0.7),
        legend.text = element_text(size=27),
        legend.title = element_markdown(size=30),
        legend.spacing.y = unit(0.1,"cm"),
        axis.title.x=element_markdown(size=35),
        axis.title.y=element_text(size=35),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_text(size=27),
        axis.text.y=element_text(size=27))
plot
 