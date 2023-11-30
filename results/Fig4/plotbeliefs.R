library(ggplot2)
library(ggtext)

library(reshape2)

#read all the files

flist<-c("ctpbeliefsL0.txt","ctpbeliefsLall.txt")

dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data<-data.frame(do.call(rbind,dt_list))
names(data)<-c("Reasoning","1","2","3","4","5")
data[1,1]<-"no **ToM**"
data[2,1]<-"**ToM**"
data <- data[c(2,1),]

#
df <- melt(data, id.vars="Reasoning")
df[,3] <- round (df[,3],2)

plot<-ggplot(data=df, aes(x=variable, y=value, fill=Reasoning)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  labs(x=expression("beliefs(*t*)"), y = "Frequency", fill="")+
  geom_text(aes(label=value), vjust=1.6, color="white",
            position = position_dodge(0.9), size=10)+
  #scale_fill_brewer(labels=c("**ToM** (*0 \u2264 k \u2264 4*)","no **ToM** (*k = 0*)"), palette="Paired", direction=-1)+
  theme_minimal()+
  scale_fill_manual(labels=c("**ToM** (*0 \u2264 k \u2264 4*)","no **ToM** (*k = 0*)"),values=c("#0072B2","#FFCC66"))+
  guides(fill = guide_legend(byrow = TRUE))+
  theme(legend.position = c(0.2, 0.85),
        legend.text = element_markdown(size=27),
        legend.title = element_text(size=27),
        legend.spacing.y = unit(1,"cm"),
        axis.title.x = element_markdown(size=35),
        axis.title.y=element_text(size=35),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_text(size=27),
        axis.text.y=element_text(size=27))
plot

