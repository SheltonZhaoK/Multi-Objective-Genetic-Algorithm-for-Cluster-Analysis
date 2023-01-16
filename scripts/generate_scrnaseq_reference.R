library(scrnabench)
dataDir = '../data'
outputDir = "../data/scrna_benchmarks_umap"
datasets = download_data(path = dataDir)
dataList = extract_datasets(datasets)
dataList = filter_data(dataList)
dataList = annotate_datasets(dataList)
dataList = run_log(dataList)
dataList = select_hvg(dataList)
dataList = scale_data(dataList)
dataList = run_pca(dataList, numComponents = 10)
dataList = run_umap(dataList, numDimensions = 10)
names = names(dataList)
for (i in 1:length(names))
{
   print(i)
   data = Seurat::Embeddings(dataList[[i]],'umap')
   file = paste(outputDir,"/" ,names[i],"_umap.csv", sep= "")
   write.csv(data,file, row.names = TRUE)
}

