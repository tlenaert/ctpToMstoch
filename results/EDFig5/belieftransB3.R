library(igraph)
library(ggplot2)
library(ggraph)
library(RColorBrewer)

nodes <- read.csv("ctpbeliefsB3.txt", header=T,  sep=",", as.is=T)
names(nodes) <- c("id","b","sd")
links <- read.csv("ctpltransitbeliefsB3.txt", header=F, sep="\t", as.is=T)
names(links) <- c("from", "to", "weight")
links[,3] <- round(links[,3],3) 
links <-links[links["weight"] >1,]


par(mar=c(0,0,0,0))
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)

pal <- brewer.pal(length(nodes$id)+2, "Greens")
V(net)$color <- pal[2:6]
V(net)$frame.color <- pal[7]
V(net)$label <- V(net)
V(net)$label.cex <- 3+ 7.5*V(net)$sd
V(net)$size <- V(net)$sd*25 + 20
V(net)$label.color <- "white"

n_bins<-10
sn_colorrange <- colorRampPalette(c("lightgray","black"))
sn_color <- sn_colorrange(n_bins)

E(net)$color<-sn_color[cut_number(E(net)$weight,n_bins)]
E(net)$width <- 1+(E(net)$weight*0.6)
E(net)$arrow.size <- 3

plot(net, layout=layout_in_circle, edge.curved=.25)

