
setwd("/Users/zhuoranx/Documents/Coursework/STAT970/final")
#################---------------simulation of swiss roll------------##############
library(snedata)
library(ggplot2)
library(plot3D)
library(pracma)
library(ggpubr)
library(gridExtra)

n <- 3000
t <- seq(1.5, 4.5,length.out=n) * pi
x <- t * cos(t)
y <- 21 * runif(n,0,1)
z <- t * sin(t)
swiss_roll <- data.frame(x = x, y = y, z = z)
eps = as.matrix(rnorm(9000,0,0.5),ncol=3)
swiss_roll = swiss_roll+eps
swiss_roll$color = t

#setwd("~/Documents/Coursework/STAT970/final")
#png("swissroll.png",width = 4500, height = 2000,res=600)
layout(matrix(c(1,2,3), ncol = 3), widths = c(1.5, 1.5, 0.2))
par(mar = c(5, 4, 4, 2) ) 
scatter3D(x = swiss_roll$x, y = swiss_roll$y, z = swiss_roll$z, 
          colvar = swiss_roll$color,
          pch = 20,  # Point shape
          cex = 0.5,  # Point size
          xlab = "X-axis", ylab = "Y-axis", zlab = "Z-axis")
mtext("a", side = 3, line = 0, at = -0.7, adj = 0, cex = 1.5)
par(mar = c(5, 4, 4, 2) +1) 
scatter2D(x = swiss_roll$x, y = swiss_roll$z,
          colvar = swiss_roll$color,
          pch = 20,  # Point shape
          cex = 0.5,  # Point size
          xlab = "X-axis", ylab = "Y-axis")
mtext("b", side = 3, line = 1, at = -25, adj = 0, cex = 1.5)
#dev.off()

#########################################swiss role in 100d  ######################################
# add zero padding to the Swiss roll 
swiss_roll_padded <- cbind(as.matrix(swiss_roll[,1:3]), matrix(0, nrow = nrow(swiss_roll), ncol = 97))
# Generate a random 10x10 orthogonal matrix
Q <- randortho(100)
#rotate the padded Swiss roll data into the 100D space using the orthogonal matrix
swiss_roll_mat <- swiss_roll_padded %*% Q
swiss_roll_mat = scale(swiss_roll_mat)
rownames(swiss_roll_mat)= paste0('samples',1:3000)
colnames(swiss_roll_mat) = paste0('features',1:100)

source("~/Documents/Coursework/STAT970/final/meta-visualization/R Codes/main_fun.R")
candidate.out = candidate.visual(swiss_roll_mat,method=c("PCA", "MDS","LEIM","Sammon","UMAP", "tSNE","PHATE"),
                                 umap.k= c(30), tsne.perplexity = c(30),phate.k = 30)

#ensemble
ensemble.out = ensemble.viz(candidate.out$embed.list, candidate.out$method_name)
ensemble.data=umap(as.dist(ensemble.out$ensemble.dist.mat),  n_neighbors = 30)
candidate.out$embed.list[[8]] = ensemble.data
candidate.out$method[8]='Ensemble'

plist = list()
for (k in 1:8){
  print(k)
  data.plot = data.frame(dim1=candidate.out$embed.list[[k]][,1], 
                         dim2=candidate.out$embed.list[[k]][,2], 
                         color=1:3000)
  p = ggplot(data.plot, aes(x=dim1, y=dim2)) + 
    geom_point(size=0.6, aes( color=color),show.legend = F) +
    ggtitle(candidate.out$method[k])+
    theme_classic()+
    scale_color_viridis_c()
    #scale_color_gradientn(colours = jet.col (n = 2500, alpha = 1))
  plist[[k]] = p
}
combined_plot <- grid.arrange(grobs = plist, ncol = 4, nrow = 2)
png("swissroll_metaviz.png",width = 4800, height = 2700,res=600)
combined_plot <- grid.arrange(grobs = plist, ncol = 4, nrow = 2)
dev.off()



#ensemble
ensemble.out = ensemble.viz(candidate.out$embed.list, candidate.out$method_name)
data.plot = data.frame(eigen.score = c(ensemble.out$eigenscore), 
                       method = rep(candidate.out$method_name, 
                                    each=dim(swiss_roll_mat)[1]))
ggplot(data.plot, aes(x=reorder(method, eigen.score, FUN=median), y=eigen.score)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("eigenscore") + xlab("method")


