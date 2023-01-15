# Tutorial: Multi-Objective Genetic Algorithm for Cluster Analysis 
This package aims to cluster single-cell RNA sequencing transcriptomes using multi-objective genetic algorithm.  
The tutorial is arranged as follows:
1. Dependency requirements
2. Data Preparation
3. Test demo

## Dependency requirements
The scripts are written in R and Python. In R, the following packages are needed: splatter, scater, scrnabench, aricode. In python, the following packages are needed: scanpy, numpy, pandas, deap, scipy, sklearn, scoop, numba, umap.umap_, validclust, seaborn, matplotlib. 

## Data Preparation
Data are generated and preprocessed in R. To generate data, please, create a directory called 'data'. If you are in "./scripts", please execute the followings.
```
mkdir ../data
cd ../data
```

To generate the datasets, please make sure your current working directory is in "data".

### scRNA-seq reference datesets
To generate scRNA-seq reference datesets for internal validation, please execute the following commands:

```
mkdir scrna_benchmarks_umap
Rscript ../scripts/generate_scrnaseq_benchmarks.R
```

### scRNA-seq synthetic datasets
To generate scRNA-seq synthetic datasets for external validation, please execute the following commands:

```
mkdir synthetic_datasets
Rscript ../scripts/generate_scrnaseq_synthetic.R
```

### scRNA-seq reference datesets under perturbations
To generate scRNA-seq reference datasets for metamorphic testing, please execute the following commands:
```
mkdir metamorphic_test
cd metamorphic_test
mkdir 1 2 3 4 5 6
Rscript ../../scripts/generate_metamorphic_datasets.R
```
