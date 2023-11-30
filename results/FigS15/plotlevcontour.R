alllibrary(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(metR)


#read all the files
flist<-list.files(pattern="ctperrorL6fine.txt")
dt_list <- lapply(flist, read.delim, header = FALSE)

#prepare data for plotting
data.cont<-data.frame(dt_list)
names(data.cont)<- c("beta", "eps", "mse","rmse","avgb","avgl","avgd","diff","avgf")
min(data.cont$rmse)
min(data.cont$mse)

#make smooth contour
mod <- loess(diff ~ beta + eps, data = data.cont,span=0.1)
df <- with(data.cont, 
           expand.grid(eps = seq(min(eps), max(eps), len = 100),
                       beta = seq(min(beta), max(beta), len = 100)))

df$yvar <- c(predict(mod, newdata = df, se = FALSE))

min(data.cont$diff)
max(data.cont$diff)


myplot <-ggplot(df, aes(x=beta, y=eps, z=yvar))+
  geom_raster(aes(fill = yvar))+
  geom_contour(color="white",bins=20)+
  scale_fill_viridis_c(limits = c(0.5,4.0), breaks = c(1.5, 2.5, 3.5)) +
  guides(fill = guide_colourbar(barheight=20,nbin=100)) +
  geom_vline(xintercept = 0.49, colour="red", linetype="dashed")+
  geom_hline(yintercept = 0.25, colour="red", linetype="dashed")+
  geom_point(x=0.49, y=0.25, size=6, colour="red")+
  annotate("text", x=0.495, y=0.26, label= "("~beta~"=0.49,"~epsilon~"=0.25)", colour="white", size=5, hjust="left")+
  geom_text_contour(aes(z = yvar),stroke = 0.2, size =6,skip=0.5)+
  xlab(beta~"(selection strength)")+
  ylab(epsilon~"(reasining error)")+
  labs(fill='Level')+ 
  #coord_fixed()+
  theme_minimal()+
  theme(#aspect.ratio = 1,
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

