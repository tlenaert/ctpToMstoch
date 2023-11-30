library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggtext)

#read all the files
flist<-list.files(pattern="ctpmbperKLall.txt")

dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data<-data.frame(dt_list)
names(data)<- c("beta", "level", "correct", "negative", "optimistic", "mb")
datasub<-select(data,"level", "mb")
df <- melt(datasub, id.vars="level")
df[,3] <- round (df[,3],2)


plot<-ggplot(data=df, aes(x=level, y=value)) +
  geom_bar(stat="identity", color="black")+
  labs(x="*k*-levels", y = "*t-T* Frequency")+
  theme_minimal()+
  ylim(-0.3,0.3)+
  guides(fill = guide_legend(byrow = TRUE))+
  theme(legend.position = "none",
        axis.title.x=element_markdown(size=35),
        axis.title.y=element_markdown(size=35),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_text(size=27),
        axis.text.y = element_text(size=27))
plot


