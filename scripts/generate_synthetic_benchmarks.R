# use R > 4.0
# to install required packages
# devtools::install_version("RcppAnnoy", "0.0.16", repos="http://cran.us.r-project.org")
# install.packages("BiocManager")
# BiocManager::install("BiocNeighbors")
# BiocManager::install("scater")
# BiocManager::install("splatter")
# devtools::install_version("RcppAnnoy")

# create synthetic data

suppressPackageStartupMessages({
  library(splatter)
  library(scater)
  library(scrnabench)
  library(aricode)
})

dirPath <- '../data/synthetic_datasets/'
set.seed(1)

numberClusters <- 1
numberCells <- 10000 #number of instances could be changed here

params <- newSplatParams()
params <- setParam(params, "batchCells", numberCells)

for (i in c(1:6))
{   
    id <- 1
    dataList <- NULL
    numberClusters = numberClusters * 2
    prob = 1/numberClusters
    
    sim <- splatSimulate(params, group.prob = rep(prob, numberClusters), method = "groups")
    dataList[[id]] <- counts(sim)
    names(dataList)[id] <- paste('synthetic_', i, sep='')
    describe(dataList)
    dataList <- filter_data(dataList)
    dataList <- run_log(dataList)
    dataList <- select_hvg(dataList)
    dataList <- scale_data(dataList)
    dataList <- run_pca(dataList, numComponents = 10)
    dataList <- run_umap(dataList, numDimensions = 10)
    dataList <- run_kmeans_clustering(dataList, reductionType = 'umap', numberClusters)
    dataList <- run_seurat_clustering(dataList, reductionType = 'umap', numberComponents = 2)

    data <- Seurat::Embeddings(dataList[[id]],'umap')
    write.csv(data, paste(dirPath, 'synthetic_embeddings_', i, '.csv', sep=''))
    cells <- rownames(data)
    labels <- as.data.frame(sim@colData)[cells, c('Batch', 'Group')]
    labels <- cbind(labels, dataList[[id]]$KMEANS_CLUSTER_UMAP, dataList[[id]]$SEURAT_CLUSTER_UMAP)
    labels$Batch <- as.numeric(as.factor(labels$Batch))
    labels$Group <- as.numeric(as.factor(labels$Group))
    labels <- cbind(labels$Batch, dataList[[id]]$SEURAT_CLUSTER_UMAP, dataList[[id]]$KMEANS_CLUSTER_UMAP, labels$Group)
    labels <- as.data.frame(labels)
    colnames(labels) <- c('did', 'seurat_cluster', 'kmeans_cluster', 'labels')
    write.table(data.frame('ID'=rownames(labels), labels), paste(dirPath, 'labels_', i, '.txt', sep=''), sep='\t', quote=F, row.names=F)
    print(run_silhouette(dataList, method='kmeans', reductionType='umap'))
    print(ARI(labels$kmeans_cluster, labels$labels))
}

