library(scrnabench)

args <- commandArgs(trailingOnly = TRUE)
seed = args[1]

datasets = load_data(path = '../data')
dataList = extract_datasets(datasets)
method = 'seurat'
transformationType = 'log'
numberClusters = 10
metamorphicTests = c(1:6)

set.seed(seed)
perturbations <- c('Permute Cells', 'Modify Gene Counts', 'Add Duplicate Cell', 'Permute Genes', 'Add Zero Variance Gene', 'Flip Gene Counts')

i <- 0
metamorphicReportList <- NULL
for(test in metamorphicTests)
{
    i <- i + 1
    metamorphicReportList[[i]] <- run_metamorphic_test(dataList, test, method, transformationType, seed, numberClusters)
}
names(metamorphicReportList) <- perturbations[metamorphicTests]

reportTable <- NULL
reportTable$Dataset = metamorphicReportList[[1]]$Dataset

i <- 1
for(test in perturbations)
{
    reportTable[test] = c(metamorphicReportList[[i]]["Silhouette UMAP"]) 
    i <- i + 1
}

print(reportTable)
path = paste("../output/seurat_metamorphic/metamorphic_log_" , seed , ".csv", sep = "")
write.csv(reportTable, file = path, row.names= F, quote = F)
