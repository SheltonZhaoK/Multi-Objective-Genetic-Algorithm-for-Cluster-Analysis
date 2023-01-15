library(scrnabench)

args <- commandArgs(trailingOnly = TRUE)

method = "seurat"

transformation = "log"

seed = args[1]

datasets = load_data(path = '../data') 

dataList = extract_datasets(datasets) 

dataList = run_clustering_workflow(dataList, method, transformation, seed = seed)

resultsTable = create_internal_cluster_validation_report(dataList, method)

path = paste("../output/seraut/",toupper(method), "_CLUSTERING_" ,toupper(transformation),seed , ".csv", sep ="") 

write.csv(resultsTable, file = path, row.names= F, quote = F)



