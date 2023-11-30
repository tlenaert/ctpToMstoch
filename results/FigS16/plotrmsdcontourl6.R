library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(metR)


#read all the files
flist<-list.files(pattern="cpieerrorLall.txt")
dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data.cont<-data.frame(dt_list)
names(data.cont)<- c("beta", "eps", "mse","rmse","avgb","avgl","avgd","diff","avgf")
min(data.cont$rmse)
min(data.cont$mse)

#make smooth contour
mod <- loess(rmse ~ beta + eps, data = data.cont, span=0.3)
df <- with(data.cont, 
           expand.grid(eps = seq(min(eps), max(eps), len = 100),
                       beta = seq(min(beta), 0.2, len = 100))) #max(beta)

df$yvar <- c(predict(mod, newdata = df, se = FALSE))

min(data.cont$rmse)
max(data.cont$rmse)

myplot <-ggplot(df, aes(x=beta, y=eps, z=yvar))+
  geom_raster(aes(fill = yvar))+
  geom_contour(color="white",bins=15)+
  scale_fill_viridis_c(limits = c(0.0,0.21), breaks = c(0.0, 0.05, 0.1, 0.15, 0.2)) +
  guides(fill = guide_colourbar(barheight=20,nbin=100)) +
  geom_vline(xintercept = 0.04, colour="red", linetype="dashed")+
  geom_hline(yintercept = 0.76, colour="red", linetype="dashed")+
  geom_point(x=0.04, y=0.76, size=6, colour="red")+
  annotate("text", x=0.045, y=0.77, label= "("~beta~"=0.04,"~epsilon~"=0.76)", colour="white", size=5, hjust="left")+
  geom_text_contour(aes(z = yvar),stroke = 0.2, size =6)+
  xlab(beta~"(selection strength)")+
  ylab(epsilon~"(reasining error)")+
  labs(fill='RMSE')+ 
#  coord_fixed()+
  theme_minimal()+
  theme(#aspect.ratio = 1,
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

