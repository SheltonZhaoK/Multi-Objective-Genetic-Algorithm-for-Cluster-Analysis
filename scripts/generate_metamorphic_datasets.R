library(scrnabench)

modify_gene_counts <- function(dataList)
{
  if(is.list(dataList)){
     for (i in (1:length(names(dataList)))) {
        geneIndex <- sample.int(nrow(dataList[[i]]), 1)
        dataList[[i]][geneIndex,] <- 2 * dataList[[i]][geneIndex,] + 1
     }}
  else
  {
    stop("A list of datasets is required to modify gene counts.")
  }
  return(dataList)
}

add_duplicate_cells <- function(dataList)
{
  if(is.list(dataList)){
     for (i in (1:length(names(dataList)))) {
        cellIndex <- sample.int(ncol(dataList[[i]]), 1)
        duplicatedCell <- dataList[[i]][, cellIndex]
        columnNames <- c(colnames(dataList[[i]]), paste(colnames(dataList[[i]])[cellIndex], '-dup', sep=''))
        dataList[[i]] <- cbind(dataList[[i]], duplicatedCell)
        colnames(dataList[[i]]) <- columnNames
     }}
  else
  {
    stop("A list of datasets is required to add duplicate cells.")
  }
  return(dataList)
}
 
add_zero_variance_gene_counts <- function(dataList)
{
  if(is.list(dataList)){
     for (i in (1:length(names(dataList)))) {
        zeroVarianceCounts <- rep(1, ncol(dataList[[i]]))
        dataList[[i]] <- rbind(dataList[[i]], zeroVarianceCounts)
     }}
  else
  {
    stop("A list of datasets is required to add zero variance gene counts.")
  }
  return(dataList)
}

flip_gene_counts <- function(dataList)
{ 
  if(is.list(dataList)){
     for (i in (1:length(names(dataList)))) {
        geneIndex <- sample.int(nrow(dataList[[i]]), 1)
        dataList[[i]][geneIndex,] <- (-1) * dataList[[i]][geneIndex,]
     }}
  else
  { 
    stop("A list of datasets is required to modify gene counts.")
  }
  return(dataList)
}

perturb_datasets <- function(dataList, perturbationType = 1)
{
    dataList <- switch(
                    perturbationType,
                    '1'= permute_columns(dataList),
                    '2'= modify_gene_counts(dataList),
                    '3'= add_duplicate_cells(dataList),
                    '4'= permute_rows(dataList),
                    '5'= add_zero_variance_gene_counts(dataList),
                    '6'= flip_gene_counts(dataList)
                    )
   return(dataList)
}

outputDir = "../data/metamorphic_test"
dataDir = '../data'
datasets = load_data(path = dataDir)

for (perturbation in 1:6)
{
   dataList = extract_datasets(datasets)
   permutedList <- perturb_datasets(dataList, perturbation)
   permutedList = filter_data(permutedList)
   permutedList = annotate_datasets(permutedList)
   permutedList = run_log(permutedList)
   permutedList = select_hvg(permutedList)
   permutedList = scale_data(permutedList)
   permutedList = run_pca(permutedList, numComponents = 10)
   permutedList = run_umap(permutedList, numDimensions = 10)
   names = names(permutedList)
   for (i in 1:length(names))
   {
      data = Seurat::Embeddings(permutedList[[i]],'umap')
      file = paste(outputDir,"/" ,perturbation,"/" ,names[i],"_umap.csv", sep= "")
      write.csv(data,file, row.names = TRUE)
   }
}
