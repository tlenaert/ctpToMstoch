library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(metR)


#read all the files
flist<-list.files(pattern="ctperrorL6.txt")
dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data<-data.frame(dt_list)
names(data)<- c("beta", "eps", "mse","rmse")
min(data$rmse)
min(data$mse)

myplot <-ggplot(data, aes(x=beta, y=eps, z=rmse))+
  geom_contour_filled(color="white",binwidth=0.025, aes(fill = stat(level)))+
  geom_vline(xintercept = 0.45, colour="white", linetype="dashed")+
  geom_hline(yintercept = 0.24, colour="white", linetype="dashed")+
  geom_point(x=0.45, y=0.24, size=6, colour="white")+
  annotate("text", x=0.455, y=0.25, label= "("~beta~"=0.45,"~epsilon~"=0.24)", colour="white", size=5, hjust="left")+
  geom_text_contour(aes(z = rmse),stroke = 0.2, size =6)+
  xlab(beta~"(selection strength)")+
  ylab(epsilon~"(reasining error)")+
  labs(fill='RMSE')+ 
  coord_fixed()+
  theme_minimal()+
  theme(aspect.ratio = 1,
        legend.text = element_text(size=18),
        legend.title = element_text(size=20),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=20),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_text(size=20),
        axis.text.y=element_text(size=20))
myplot

