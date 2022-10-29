library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)

#read all the files
flist<-list.files(pattern="ctpstepsexpl4.txt")
dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data<-data.frame(dt_list)
names(data)<- c("beta", "eps","S0", "S1", "S2", "S3", "S4")
data<-data[,-1]
data<-data[,-1]
mkpavg <- c(0.071, 0.356, 0.370, 0.153, 0.049)


data<-rbind(data, mkpavg)

newdata <-c(1,2,3,4,5)
newdata <-rbind(newdata,data)

newdata<-t(newdata)
newdata<-data.frame(newdata)
names(newdata)<-c("steps", "Present model", "Experiment")
df <- melt(newdata, id.vars="steps")

plot<-ggplot(data=df, aes(x=steps, y=value, fill=variable)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  labs(x="Steps", y = "Frequency", fill="")+
  scale_fill_manual(values=c("#0072B2","#FF6600"))+
  scale_x_continuous(breaks=0:6)+
  theme_minimal()+
  expand_limits(y=0.4)+
  theme(legend.position = c(0.8, 0.8),
        legend.text = element_text(size=27),
        legend.title = element_text(size=27),
        axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_text(size=27),
        axis.text.y=element_text(size=27))
plot

