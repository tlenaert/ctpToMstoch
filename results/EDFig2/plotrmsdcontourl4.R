library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(metR)


#read all the files
flist<-list.files(pattern="ctperrorL4.txt")
dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data<-data.frame(dt_list)
names(data)<- c("beta", "eps", "mse","rmse","avgb","avgl","avgd","diff","avgf")
min(data$rmse)
min(data$mse)

myplot <-ggplot(data, aes(x=beta, y=eps, z=rmse))+
  geom_contour_filled(color="white",binwidth=0.05, aes(fill = stat(level)))+
  geom_vline(xintercept = 0.3, colour="white", linetype="dashed")+
  geom_hline(yintercept = 0.18, colour="white", linetype="dashed")+
  geom_point(x=0.3, y=0.18, size=6, colour="white")+
  annotate("text", x=0.305, y=0.19, label= "("~beta~"=0.3,"~epsilon~"=0.18)", colour="white", size=7, hjust="left")+
  geom_text_contour(aes(z = rmse),stroke = 0.2, size =7)+
  xlab(beta~"(selection strength)")+
  ylab(epsilon~"(reasoning error probability)")+
  labs(fill='Average level')+ 
  #coord_fixed()+
  theme_minimal()+
  theme(
        legend.position = "none",
        legend.text = element_text(size=18),
        legend.title = element_text(size=27),
        axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_text(size=27),
        axis.text.y=element_text(size=27))
myplot


