#####simulating a tetrahedron
library(ggplot2)
library(plot3D)
library(pracma)
library(ggpubr)
library(gridExtra)

set.seed(42)
generate_tetra = function(n = 500, sd = 0.05){
  A<-c(0,0,0)
  B<-c(1,0,0)
  C<-c(1/2,sqrt(3)/2,0)
  D<-c(1/2,sqrt(3)/6,sqrt(6)/3)
  cluster_A = matrix(rnorm(3*n,sd=sd),ncol = 3)+matrix(rep(matrix(A,ncol=3),times=n),ncol=3,byrow = T)
  cluster_B = matrix(rnorm(3*n,sd=sd),ncol = 3)+matrix(rep(matrix(B,ncol=3),times=n),ncol=3,byrow = T)
  cluster_C = matrix(rnorm(3*n,sd=sd),ncol = 3)+matrix(rep(matrix(C,ncol=3),times=n),ncol=3,byrow = T)
  cluster_D = matrix(rnorm(3*n,sd=sd),ncol = 3)+matrix(rep(matrix(D,ncol=3),times=n),ncol=3,byrow = T)
  X = rbind(cluster_A,cluster_B,cluster_C,cluster_D)
  return(X)
}

X = generate_tetra(500,sd=0.1)
png("tetrahedron.png",width = 3000, height = 3000,res=600)
scatter3D(x = X[,1], y = X[,2], z = X[,3], 
          colvar = rep(1:4,each=500),
          col = c("#d73027", "#fee090", "#abdda4", "#4575b4"), 
          pch = 20,  # Point shape
          cex = 0.5,  # Point size
          xlab = "X-axis", ylab = "Y-axis", zlab = "Z-axis",
          clab=NULL,colkey=F)
rownames(X)= paste0('samples',1:2000)
colnames(X) = paste0('features',1:3)
dev.off()


source("~/Documents/Coursework/STAT970/final/meta-visualization/R Codes/main_fun.R")
candidate.out = candidate.visual(X,method=c("PCA", "MDS","LEIM","Sammon","UMAP", "tSNE","PHATE"),
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
                         color=as.factor(rep(1:4,each=500)))
  p = ggplot(data.plot, aes(x=dim1, y=dim2)) + 
    geom_point(size=0.6, aes( color=color),show.legend = F) +
    ggtitle(candidate.out$method[k])+
    theme_classic()+
    scale_color_manual(values = c("#d73027", "#fee090", "#abdda4", "#4575b4"))
  plist[[k]] = p
}
png("tetrahedron_viz.png",width = 4800, height = 2700,res=600)
combined_plot <- grid.arrange(grobs = plist, ncol = 4, nrow = 2)
dev.off()



#ensemble---eigenscore
data.plot = data.frame(eigen.score = c(ensemble.out$eigenscore), 
                       method = rep(candidate.out$method_name,each=dim(X)[1]))
ggplot(data.plot, aes(x=reorder(method, eigen.score, FUN=median), y=eigen.score)) +
  theme_classic()+
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("eigenscore") + xlab("method")
ggsave('tetrahedron_eigenscore.png')

