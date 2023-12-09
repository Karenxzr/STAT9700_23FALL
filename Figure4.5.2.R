library(rARPACK)
library(MASS)
library(dimRed)
library(uwot)
library(cluster)
library(phateR)
library(Rtsne)
library(stringr)
library(ggplot2)
library(dplyr)
library(gridExtra)
##############################discrete data##########################################
load('sim.discrete.RData')
data_discrete_pc = read.csv('pca_embeddings_discrete_pbmc3k.csv')
rownames(data_discrete_pc) = paste0('cell',1:nrow(data_discrete_pc))
source("~/Documents/Coursework/STAT970/final/meta-visualization/R Codes/main_fun.R")
candidate.out = candidate.visual(data_discrete_pc[,1:20],
                                 method=c("PCA", "MDS","LEIM","Sammon","UMAP", "tSNE","PHATE"),
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
                         color=as.factor(unlist(lapply(str_split(rownames(sim.discrete@meta.data),'_'),
                                                       function(x)x[[2]]))))
  p = ggplot(data.plot, aes(x=dim1, y=dim2)) + 
    geom_point(size=0.6, aes( color=color),show.legend = F) +
    ggtitle(candidate.out$method[k])+
    theme_classic()+
    scale_color_manual(values = c( "#006d2c","#d73027","#8856a7","#a1d99b","#4575b4"))
  plist[[k]] = p
}

png("discrete_cluster_viz.png",width = 4800, height = 2700,res=600)
combined_plot <- grid.arrange(grobs = plist, ncol = 4, nrow = 2)
dev.off()

#ensemble
data.plot = data.frame(eigen.score = c(ensemble.out$eigenscore), 
                       method = rep(candidate.out$method_name, 
                                    each=dim(data_discrete_pc)[1]))
ggplot(data.plot, aes(x=reorder(method, eigen.score, FUN=median), y=eigen.score)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("eigenscore") + xlab("method")
ggsave('discrete_eigenscore.png')


##############################continuous data##########################################
setwd("~/Documents/Coursework/STAT970/final/data")
load('sim.continuous.RData')
data_continuous_pc = read.csv('pca_embeddings_continuous_paul.csv')
rownames(data_continuous_pc) = paste0('cell',1:nrow(data_continuous_pc))
source("~/Documents/Coursework/STAT970/final/meta-visualization/R Codes/main_fun.R")

candidate.out = candidate.visual(data_continuous_pc,
                                 method=c("PCA", "MDS","LEIM","Sammon","UMAP", "tSNE","PHATE"),
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
                         color=as.factor(sim.continuous$group))
  p = ggplot(data.plot, aes(x=dim1, y=dim2)) + 
    geom_point(size=0.6, aes( color=color),show.legend = F) +
    ggtitle(candidate.out$method[k])+
    theme_classic()+
    scale_color_manual(values = c("#d73027", "#abdda4","#66c2a4", "#006d2c"))
  plist[[k]] = p
}
png("trajectory_viz.png",width = 4800, height = 2700,res=600)
combined_plot <- grid.arrange(grobs = plist, ncol = 4, nrow = 2)
dev.off()

#ensemble
data.plot = data.frame(eigen.score = c(ensemble.out$eigenscore), 
                       method = rep(candidate.out$method_name, 
                                    each=dim(data_continuous_pc)[1]))
ggplot(data.plot, aes(x=reorder(method, eigen.score, FUN=median), y=eigen.score)) +
  geom_boxplot(outlier.size = 0.5) + theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1)) +
  ylab("eigenscore") + xlab("method")
ggsave('trajectory_eigenscore.png')
