
library(ggplot2)
library(igraph)
library(ggraph)
library(RColorBrewer)

nodes <- read.csv("ctpnodesB031.txt", header=T, as.is=T)
links <- read.csv("ctpedgesB031.txt", header=T, as.is=T)
links[,3] <- round(links[,3],3) 
links <-links[links["weight"] >8,]

par(mar=c(0,0,0,0))
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
strength(net,mode="in")
strength(net,mode="out")

group_color_nodes <-  c("#D55E00","#E69F00","#009E73","#0072B2","#000000")
V(net)$color <- adjustcolor(group_color_nodes[V(net)$k+1],1.0)
V(net)[24]$color <- adjustcolor(group_color_nodes[V(net)[24]$k+1])
V(net)$frame.color <- "#000000"
V(net)$frame <- 5
V(net)$label <- V(net)$b+1
V(net)$label.cex <- 2+ 15*V(net)$sd
V(net)$size <- V(net)$sd*50 + 15
V(net)$label.color <- "white"


n_bins<-10
sn_colorrange <- colorRampPalette(c("lightgray","black"))
sn_color <- sn_colorrange(n_bins)

E(net)$color<-sn_color[cut_number(E(net)$weight,n_bins)]
E(net)$width <-5
E(net)$arrow.size <- 2


l1 <- layout_on_grid(net)
l2 <- layout_with_graphopt(net, charge=1.0)
l3<- layout_on_sphere(net)
plot(net, layout=l1, edge.curved=.2)

legend(x=-1.7, y=1.0, c("k=0","k=1", "k=2", "k=3", "k=4"), pch=21,
       col="#777777", pt.bg=adjustcolor(group_color_nodes, alpha=1.0), pt.cex=7, cex=3.0, bty="n", ncol=1)


