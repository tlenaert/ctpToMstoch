library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(metR)


#read all the files
flist<-list.files(pattern="ctperrorL4fine.txt")
dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data.cont<-data.frame(dt_list)
names(data.cont)<- c("beta", "eps", "mse","rmse","avgb","avgl","avgd","diff","avgf")
min(data.cont$rmse)
min(data.cont$mse)
#make smooth contour
mod <- loess(rmse ~ beta + eps, data = data.cont, span=0.2)
df <- with(data.cont, 
           expand.grid(eps = seq(min(eps), max(eps), len = 100),
                       beta = seq(min(beta), max(beta), len = 100)))

df$yvar <- c(predict(mod, newdata = df, se = FALSE))

min(data.cont$rmse)
max(data.cont$rmse)

myplot <-ggplot(df, aes(x=beta, y=eps, z=yvar))+
  geom_raster(aes(fill = yvar))+
  geom_contour(color="white",bins=20)+
  scale_fill_viridis_c(limits = c(0.0,0.45), breaks = c(0.0, 0.1, 0.2, 0.3, 0.4)) +
  guides(fill = guide_colourbar(barheight=20,nbin=100)) +
  geom_vline(xintercept = 0.31, colour="red", linetype="dashed")+
  geom_hline(yintercept = 0.19, colour="red", linetype="dashed")+
  geom_point(x=0.31, y=0.19, size=6, colour="red")+
  annotate("text", x=0.315, y=0.2, label= "("~beta~"=0.31,"~epsilon~"=0.19)", colour="white", size=7, hjust="left")+
  geom_text_contour(aes(z = yvar),stroke = 0.2, size =7)+
  xlab(beta~"(selection strength)")+
  ylab(epsilon~"(reasoning error probability)")+
  labs(fill='Average level')+ 
  #coord_fixed()+
  theme_minimal()+
  theme(
        #legend.position = "none",
        legend.text = element_text(size=18),
        legend.title = element_text(size=27),
        axis.title.x=element_text(size=35),
        axis.title.y=element_text(size=35),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_text(size=27),
        axis.text.y=element_text(size=27))
myplot


