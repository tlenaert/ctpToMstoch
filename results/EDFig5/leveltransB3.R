library(igraph)
library(ggplot2)
library(ggraph)
library(RColorBrewer)

nodes <- read.csv("ctplevelsB3.txt", header=T,  sep=",", as.is=T)
names(nodes) <- c("id","k","sd")
links <- read.csv("ctptransitlevelsB3.txt", header=F, sep="\t", as.is=T)
names(links) <- c("from", "to", "weight")
links[,3] <- round(links[,3],1) 
links <-links[links["weight"] > 1,]


par(mar=c(0,0,0,0))
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
group_color_nodes <-  c("#D55E00","#E69F00","#009E73","#0072B2","#000000")
V(net)$color <- group_color_nodes[V(net)]
V(net)$frame.color <- group_color_nodes[V(net)]
V(net)$label <- nodes$k
V(net)$label.cex <- 3+ 7.5*V(net)$sd
V(net)$size <- V(net)$sd*25 + 20
V(net)$label.color <- "white"

n_bins<-10
sn_colorrange <- colorRampPalette(c("lightgray","black"))
sn_color <- sn_colorrange(n_bins)

E(net)$color<-sn_color[cut_number(E(net)$weight,n_bins)]
E(net)$width <- 2+(E(net)$weight*0.2)
E(net)$arrow.size <- 3

plot(net, layout=layout_in_circle, edge.curved=.25)

