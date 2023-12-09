

#######################################################################################################
###############simulation of discrete clustered single cell data based on pbmc3k dataset###############
#######################################################################################################

#### Generate discrete clusters using Splatter ####
library(splatter)
library(Seurat)
library(SeuratObject)
library(SeuratData)
library(scRNAseq)
data("pbmc3k")
pbmc3k = UpdateSeuratObject(pbmc3k)
## Fit params from pbmc3k dataset
splat.seed <- 3250879
init.params <- newSplatParams(group.prob = c(0.45, 0.15, 0.15, 0.15, 0.1),
                              de.prob = c(0.3, 0.15, 0.15, 0.15, 0.3),
                              de.facLoc = 0.3, de.facScale = 0.6,
                              seed = splat.seed)
est.params <- splatEstimate(as.matrix(pbmc3k[['RNA']]@counts), params = init.params)
## Generate synthetic data using fitted parameters
groups <- splatSimulateGroups(est.params, verbose = F)
counts <- as(groups@assays@data$counts, "matrix")
metadata <- groups@colData@listData
colnames(counts) <- paste(colnames(counts), metadata$Group, sep = "_")
#--------preprocessing---------------#
sim.discrete = SeuratObject::CreateSeuratObject(counts)
sim.discrete <- NormalizeData(sim.discrete)
sim.discrete <- FindVariableFeatures(sim.discrete, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sim.discrete)
sim.discrete <- ScaleData(sim.discrete, features = all.genes)
#get the normalized data matrix
sim.discrete.scale = sim.discrete@assays$RNA$scale.data
####################
setwd("~/Documents/Coursework/STAT970/final/data")
save(sim.discrete, file = 'sim.discrete.RData')
sim.discrete <- RunPCA(sim.discrete, features = VariableFeatures(object = sim.discrete))
sim.discrete <- FindNeighbors(sim.discrete,reduction = "pca")
pca_embeddings = Embeddings(sim.discrete,reduction = "pca")
write.csv(pca_embeddings,"pca_embeddings_discrete_pbmc3k.csv",row.names = F)
sim.discrete <- RunUMAP(sim.discrete,reduction = "pca",dims = 1:20)
DimPlot(sim.discrete, reduction = "umap")
#######################################################################################################
###############.    simulation of continuous cell cycle data based on hemato dataset  ###############
#######################################################################################################
library(Seurat)
#setwd("/Users/zhuoranx/Documents/Coursework/STAT970/final/data/cell_cycle_vignette_files")
#exp.mat <- read.table(file = "nestorawa_forcellcycle_expressionMatrix.txt",
#                      header = TRUE, as.is = TRUE, row.names = 1)
#index = sample(1:nrow(exp.mat),0.2*nrow(exp.mat))

## Load the PAUL HEMOTO dataset
paul <- PaulHSCData(ensembl = T)
index = colData(paul)$Seq_batch_ID=="SB19"
counts <- as.matrix(assay(paul, "counts"))[,index]
splat.seed <- 350879
init.params <- newSplatParams(seed = splat.seed,
                              group.prob = c(0.3,0.3,0.2,0.2),
                              de.prob = c(0.4,0.2,0.2,0.2), de.facLoc = 0.4, de.facScale = 0.6,
                              path.from = c(0,1,1,1), 
                              path.length = c(300,200,200,200),
                              path.skew = c(0.5,0.5,0.5,0.5))
splat.params <- splatEstimate(counts, params = init.params)
## Generate synthetic data using fitted parameters
groups <- splatSimulateGroups(splat.params, verbose = F)
counts <- as(groups@assays@data$counts, "matrix")
metadata <- groups@colData@listData
colnames(counts) <- paste(colnames(counts), metadata$Group, sep = "_")
#--------preprocessing---------------#
sim.continuous = SeuratObject::CreateSeuratObject(counts)
sim.continuous@meta.data$group = colData(groups)$Group
sim.continuous <- NormalizeData(sim.continuous)
sim.continuous <- FindVariableFeatures(sim.continuous, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(sim.continuous)
sim.continuous <- ScaleData(sim.continuous, features = all.genes)
####################
setwd("~/Documents/Coursework/STAT970/final/data")
save(sim.continuous, file = 'sim.continuous.RData')
sim.continuous <- RunPCA(sim.continuous, features = VariableFeatures(object = sim.continuous))
sim.continuous <- FindNeighbors(sim.continuous,reduction = "pca")
pca_embeddings = Embeddings(sim.continuous,reduction = "pca")
write.csv(pca_embeddings,"pca_embeddings_continuous_paul.csv",row.names = F)


