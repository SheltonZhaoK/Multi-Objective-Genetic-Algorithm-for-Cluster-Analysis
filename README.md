# Tutorial: Multi-Objective Genetic Algorithm for Cluster Analysis 
This package aims to cluster single-cell RNA sequencing transcriptomes using multi-objective genetic algorithm.

The tutorial is arranged as follows:
1. Dependency requirements
2. Data Preparation
3. Test demo

## Dependency requirements
The scripts are written in R and Python. In R, the following packages are needed: matrix, splatter, scater, scrnabench, and aricode. In python, the following packages are needed: scanpy, numpy, pandas, deap, scipy, sklearn, scoop, numba, umap.umap_, validclust, seaborn, and matplotlib. 

## Data Preparation
Data are generated and preprocessed in R. To generate data, please, create a directory called 'data'. If you are in "./scripts", please execute the followings.
```
mkdir ../data
cd ../data
```

To generate the datasets, please make sure your initial current working directory is always in "data".

### scRNA-seq reference datesets
To generate scRNA-seq reference datesets for internal validation, please execute the following commands in "./data":

```
mkdir scrna_benchmarks_umap
Rscript ../scripts/generate_scrnaseq_reference.R
```

### scRNA-seq synthetic datasets
To generate scRNA-seq synthetic datasets for external validation, please execute the following commands in "./data":

```
mkdir synthetic_datasets
Rscript ../scripts/generate_scrnaseq_synthetic.R
```

### scRNA-seq reference datesets under perturbations
To generate scRNA-seq reference datasets for metamorphic testing, please execute the following commands in "./data":
```
mkdir metamorphic_test
cd metamorphic_test
mkdir 1 2 3 4 5 6
cd ..
Rscript ../scripts/generate_metamorphic_datasets.R
```

## Test Demo
Three sample executions files are provided in "./scripts".

All experiments related to Seurat clustering are implemented in R, and other algorithms are implemented in Python. If you are interested how the algorithms are implemented, please look at "./scripts/clusterer.py"

To run Seurat clustering for internal validation. Please execute the following command in "./scripts". The labels computed by Seurat for external validation are generated when creating synthetic datasets. Please check the labels file.
```
Rscript test_clustering_seurat.R
```

Metamorphic testing in Seurat uses a different workflow. Please execute the following command in "./scripts".
```
Rscript test_metamorphic_seurat.R
```

To run MOGA, SOGA, Kmeans, and Scanpy for clustering. Please execute the following command in "./scripts". The example presented is external validation, and user can change the code with different input data to perform external validation and metamorphic testing.
```
python3 test_scrnaseq_benchmarks.py
```