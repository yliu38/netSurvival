# netSurvival

## Table of Contents
- **[Introduction](#introduction)**
- **[Preparations](#preparations)**
- **[Preprocessing](#preprocessing)**
- **[Data Split and Normalizations](#data-split-and-normalizations)**
- **[Univariate Log-rank Test](#univariate-log-rank-test)**
- **[One-dimensional Hierarchical Clustering](#one-dimensional-hierarchical-clustering)**
- **[Neo4j Graph Network Building and Graph Algorithm Implementation](#neo4j-graph-network-building-and-graph-algorithm-implementation)**
- **[Feature Selection](#feature-selection)**


## Introduction
netSurvival is a novel computational framework that integrates network-based feature selection with survival analysis to identify robust biomarkers from high-dimensional genomic data. Developed specifically to address the challenges of prognostic modeling in cancer genomics, netSurvival leverages the power of gene-gene interaction networks and random walk algorithms to discover biologically meaningful markers associated with patient outcomes.

<p align="center">
<img src="https://github.com/yliu38/netSurvival/blob/main/trajectory.png" width="750">
</p>

## Preparations
### R packages installation

**R code:**
``` r
# List of required packages
packages <- c(
  "dendextend", "tidyverse", "tibble", "gsubfn", "venn",
  "readxl", "data.tree", "survival", "purrr", "reticulate")

# Install missing packages
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(packages[!installed])
}

# Load the packages
lapply(packages, library, character.only = TRUE)
```

### Python packages installation

**Bash code:**
```bash
# Python 3.8 or higher
pip install py2neo numpy pandas scikit-survival lifelines
```

### Neo4j desktop setup
Please google neo4j desktop to download the neo4j software or you can use institutional neo4j server or neo4j clound

Below is the step-by-step instructions for using Neo4j desktop:

Open the neo4j software --> click "new" --> Create project --> Add Local DBMS, input password and create --> click the project made and install Plugins of APOC and Graph Data Science Library


## Preprocessing
<img src="https://github.com/yliu38/EGNF/blob/main/image/example_expression_matrix.png" width="380">

**R code:**
``` r
# load libraries
library(dendextend)
library(tidyverse)
library(tibble)
library(gsubfn)
library(readxl)
library(data.tree)
library(survival)

library(purrr)
library(reticulate)


source("https://github.com/yliu38/netSurvival/blob/main/R/functions.R")
# remove genes with 80% zeroes and na rows
exp <- remove_sparse_rows(exp)
```

## Data Split and Normalizations

**R code:**
``` r
set.seed(123)
# number of sample IDs
sample_names <- colnames(exp)
n_spl = dim(exp)[2]

# Randomly assign each patient to one of 5 folds
fold_assignments <- sample(rep(1:5, length.out = n_spl))

# Create a list to store each train/test split
cv_folds <- list()

# For each fold
for (i in 1:5) {
  # Test set: samples in the current fold
  test_ind <- sample_names[fold_assignments == i]
  
  # Training set: all other patients
  train_ind <- sample_names[fold_assignments != i]
  
  cv_folds[[i]] <- list(train = train_ind, test = test_ind)
}

# Example: Get expression matrices for fold 5, you need to run for different folds
train_ind <- cv_folds[[5]]$train
exp_train <- exp[, train_ind]
exp_test <- exp[, !colnames(exp) %in% train_ind]

# log2 and z-score normalization
# nor has options including "two.end", "up", "down" for choosing both high and low or high only or low only expressed clusters
exp_train <- norm_dat(exp_train, nor="two.end")
exp_test <- norm_dat(exp_test, nor="two.end")
```
If you have sample replicates or paired-samples, you may want to select patients for data split to avoid potential data leakage.

## Univariate Log-rank Test
**R code:**
``` r
## add blank columns (NA values) using mutate and across
## clinical_sub is a patient clinical data including case_overall_survival_mo, censoring and samplename
pt_train <- clinical_sub[clinical_sub$sample_full %in% train_ind,]
ori_ncol <- ncol(pt_train)
# reformat tables because of duplicated pid
pt_train <- pt_train %>%
  mutate(!!!setNames(rep(list(NA_real_), nrow(exp_train)), rownames(exp_train)))

for (i in seq(nrow(pt_train))) {
  tmp <- match(pt_train[i,"samplename"],colnames(exp_train))
  pt_train[i,(ori_ncol+1):ncol(pt_train)] <- exp_train[,tmp]
}

mids <- sapply((ori_ncol+1):ncol(pt_train), function(x) {
  median(pt_train[,x])
})

## replace rec values with group info
for (j in (ori_ncol+1):ncol(pt_train)) {
  pt_train[,j] <- ifelse(pt_train[,j]>mids[j-ori_ncol],"high","low")
}

## log2-rank test
pt_train$case_overall_survival_mo <- as.numeric(pt_train$case_overall_survival_mo)
pt_train$censoring <- as.numeric(pt_train$censoring)

p_values <- sapply((ori_ncol+1):ncol(pt_train), function(x) {survdiff(Surv(case_overall_survival_mo, censoring) ~ get(colnames(pt_train)[x]), data=pt_train)$pvalue})
adj_p.values <- p.adjust(p_values, method = "bonferroni")
length(adj_p.values[adj_p.values<0.05])

markers <- data.frame(
  col_name <- colnames(pt_train)[(ori_ncol+1):ncol(pt_train)],
  adj_p <- adj_p.values
)
colnames(markers) <- c("gene", "adj.p")
markers <- markers[order(markers$adj.p, decreasing = F),] 
marker <- markers$gene[1:1000]
```

## One-dimensional Hierarchical Clustering
### Output csv files for network construction

**R code:**
``` r
# directory is the location storing results, example can be "./folder_name/train_gene_"
make_tree(exp_train[marker,], directory)

# generate url file for generating nodes in Neo4j
gene_names <- rownames(exp_train)
file1 <- paste0(directory,gene_names,".csv") 
url = c("URL", file1)
write.table(url,"url_train.csv", sep=",",  col.names=F, row.names = F)
```
**please move the generated data_train folder and url_train.csv to the import folder of Neo4j**

## Neo4j Graph Network Building and Graph Algorithm Implementation
Open the neo4j software --> click the project made --> click the "..." on the right --> Open floder Import --> move the files including url_train.csv, folder for hierarchical trees to the import directory

Open terminal, run python scripts
### Build networks and implement graph-based algorithms

**Bash code:**
```python
python create_filenodes.py # creating nodes for making graph nodes
python create_nodes_sur.py # making nodes and delete file nodes
python create_relationships.py # making edges between expression events, the default cutoff for common samples shared across edges is 4, may change according to your sample size
python make_pt_nodes_rec.py # making nodes of patients
python make_relationships_pt_exp.py # making edges between patient nodes and expression event nodes
python output_id_table.py # output node ids for following feature selection process (id_gene_map.csv)

# after database construction, run random walk algorithm
python ran_walk.py
```
## Feature Selection
Please run DB_random_walk_results.R, change input and output based on your real situations.

