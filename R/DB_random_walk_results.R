rm(list=ls())
gc()
options(stringsAsFactors = F)

# Load library
library(dplyr)
library(purrr)
library(reticulate)
library(venn)

# load data
map_dat <- read.csv("id_gene_map.csv")
pd <- import("pandas")

# read files
all_files <- list.files(path = "../algorithm_res/",pattern = "[0-9]+\\.pickle")
file_sig <- paste0("../algorithm_res/",all_files)

var_names <- paste0("pickle_",sapply(strsplit(all_files,"_|\\."),function(x)x[2]))

# Assign values to variables
for (i in seq_along(var_names)) {
  assign(var_names[i], pd$read_pickle(file_sig[i]))
}

# control type1 error
all_p <- c()
for(i in seq(length(var_names))) {
  all_p <- c(all_p, unlist(get(var_names[i])[[1]][1]))
}

all_p.adj <- p.adjust(all_p, method = "fdr")
all_p.adj <- all_p.adj[all_p<0.05]
table(all_p<0.05); table(all_p.adj<0.05)

# combine the results
sig_paths <- c()
sig_kms <- c()
sig_kms_not_in_path <- c()
patients_in_each_path <- c()

for(i in seq(length(var_names))) {
  sig_paths <- c(sig_paths, get(var_names[i])[[1]][[2]][[1]])
  sig_kms <- c(sig_kms, get(var_names[i])[[1]][[2]][[2]])
  sig_kms_not_in_path <- c(sig_kms_not_in_path, get(var_names[i])[[1]][[2]][[3]])
  patients_in_each_path <- c(patients_in_each_path, get(var_names[i])[[1]][[2]][[4]])
}

# check maximum number of pts in paths
pt_size <- unlist(lapply(patients_in_each_path, length))
max(pt_size)

# remove duplicated lists (almost none)
rm_index <- which(duplicated(sig_paths))

# filtering paths by pt size and adjusted p-value
index1 <- which(all_p.adj < 0.05)
index2 <- which(pt_size >= 10)
index_f <- intersect(index1, index2)
index_f <- setdiff(index_f,rm_index)

sig_paths <- sig_paths[index_f]
sig_kms <- sig_kms[index_f]
sig_kms_not_in_path <- sig_kms_not_in_path[index_f]
patients_in_each_path <- patients_in_each_path[index_f]

pvalues <- all_p.adj[index_f]
pt_size <- pt_size[index_f]

# add column for better survival / worse survival
## function to calculate area under the curve using Trapezoid Rule
cal_auc <- function(x, y) {
  auc = 0
  for(i in 2:nrow(x))
    h = x[i,1] - x[i-1,1]
    auc = auc + h * (y[i-1,1] + y[i,1]) / 2
    return(auc)
}

survival_comparison <- c()
for(i in seq(length(sig_paths))){
  auc1 = cal_auc(sig_kms[[i]]['duration'], sig_kms[[i]]['survival_probability'])
  auc2 = cal_auc(sig_kms_not_in_path[[i]]['duration'],sig_kms_not_in_path[[i]]['survival_probability'])
  survival_comparison <- c(survival_comparison, ifelse(auc1<auc2,"worse","better"))
}

# output path results
df_path <- data.frame(
  path = unlist(lapply(sig_paths,function(x) paste(x,collapse = ", "))),
  pvalue = pvalues,
  patients = unlist(lapply(patients_in_each_path,function(x) paste(x,collapse = ", "))),
  pt_size = pt_size,
  survival_comparison = survival_comparison
)
colnames(df_path) <- c("path","pvalue","patients", "pt_size", "survival_comparison")
df_path$length <- sapply(sig_paths, function(x)length(x))
## convert node to genes
df_path$gene <- ""
for(i in seq(nrow(df_path))) {
  tmp <- unlist(sig_paths[[i]])
  tmp2 <- map_dat[match(tmp,map_dat$nodeId),"name"]
  df_path$gene[i] <- paste(tmp2,collapse = ", ")
}
table(df_path$survival_comparison)
## order
df_path <- df_path[order(df_path$pvalue, decreasing = F),]

write.csv(df_path, "../survival/WT_fre_sig_path_042325.csv")

# do not separate
# remove duplicates in each path and count events
for(i in seq(length(sig_paths))) {
  tmp <- unlist(sig_paths[[i]])
  sig_paths[[i]] <- unique(tmp)
}
combined_all <- unlist(flatten(sig_paths))
sig_fre_all <- as.data.frame(table(combined_all))


# do the same for all paths
num_index <- sapply(strsplit(all_files,"_|\\."),function(x)x[2])
file_allpath <- paste0("../algorithm_res/test_",num_index,"_allpath.pickle")
var_names <- paste0("pickle_",num_index)

## Assign values to variables
for (i in seq_along(file_allpath)) {
  assign(var_names[i], pd$read_pickle(file_allpath[i]))
}

allpath <- c()
for(i in seq(length(var_names))) {
  allpath <- c(allpath, get(var_names[i]))
}
## remove duplicated lists (almost none)
allpath_unique <- unique(allpath)
length(allpath)
length(all_p)

## count
allpath_unique_uni <- allpath_unique
for(i in seq(nrow(allpath_unique_uni))) {
  tmp <- unlist(allpath_unique_uni[[i]])
  allpath_unique_uni[[i]] <- unique(tmp)
}
combined <- unlist(flatten(allpath_unique_uni))
all_fre <- as.data.frame(table(combined))


sig_fre_all$gene <- map_dat[match(sig_fre_all[,1],map_dat$nodeId),"name"]
all_fre$gene <- map_dat[match(all_fre[,1],map_dat$nodeId),"name"]

# aggregate
all_fre_new <- aggregate(all_fre$Freq, by=list(gene=all_fre$gene), FUN=sum)
## all
sig_fre_new_a <- aggregate(sig_fre_all$Freq, by=list(gene=sig_fre_all$gene), FUN=sum)
sig_fre_new_a$all <- all_fre_new[match(sig_fre_new_a$gene,all_fre_new$gene),2]


# fisher exact test 
num_allpath_unique <- length(allpath_unique)
## all
loc <- length(sig_paths)
sig_fre_new_a$p_value <- ""
for(i in seq(nrow(sig_fre_new_a))){
  dat <- matrix(c(sig_fre_new_a[i,2],loc-sig_fre_new_a[i,2],sig_fre_new_a[i,3]-sig_fre_new_a[i,2],num_allpath_unique+loc-sig_fre_new_a[i,2]-loc-sig_fre_new_a[i,3]-loc), nrow=2)
  sig_fre_new_a$p_value[i] <- fisher.test(dat, alternative = "greater")$p.value
}

sig_fre_new_a$adj.pvalue <- p.adjust(sig_fre_new_a$p_value, method = "fdr")
colnames(sig_fre_new_a)[2:3] <- c("sig_fre", "all_fre") 
sig_fre_new_a <- sig_fre_new_a[order(sig_fre_new_a$adj.pvalue, decreasing = F),]
length(which(sig_fre_new_a$adj.pvalue<0.05))
write.csv(sig_fre_new_a, "../survival/sig_fre_all_042325.csv")


# check overlaps over 5 folds
sig_gene1 <- read.csv("../survival/sig_fre_all_031725.csv")
sig_gene2 <- read.csv("../survival/sig_fre_all_041925.csv")
sig_gene3 <- read.csv("../survival/sig_fre_all_042025.csv")
sig_gene4 <- read.csv("../survival/sig_fre_all_042225.csv")
sig_gene5 <- read.csv("../survival/sig_fre_all_042325.csv")


# Extract significant genes
g1 <- sig_gene1 %>% filter(adj.pvalue < 0.05) %>% pull(gene)
g2 <- sig_gene2 %>% filter(adj.pvalue < 0.05) %>% pull(gene)
g3 <- sig_gene3 %>% filter(adj.pvalue < 0.05) %>% pull(gene)
g4 <- sig_gene4 %>% filter(adj.pvalue < 0.05) %>% pull(gene)
g5 <- sig_gene5 %>% filter(adj.pvalue < 0.05) %>% pull(gene)

# Create named list
venn_input <- list(
  Fold1 = g1,
  Fold2 = g2,
  Fold3 = g3,
  Fold4 = g4,
  Fold5 = g5
)

# Plot Venn with beautified appearance
venn_plot <- venn(
  venn_input,
  ilabels = "counts",
  zcolor = "style",
  ilcs = 1.2,
  sncs = 1.2,
  snames = c("Fold 1", "Fold 2", "Fold 3", "Fold 4", "Fold 5"),
  box = FALSE,
  ggplot = TRUE
)

# Save to file (high-res PNG)
ggplot2::ggsave("../survival/venn_5fold_genes_adj_pval.png", venn_plot, width = 8, height = 6, dpi = 300)

# Count genes occurring in at least 3 folds
all_genes <- c(g1, g2, g3, g4, g5)
gene_freq <- table(all_genes)
selected_genes <- names(gene_freq[gene_freq >= 3])

# Output dataframe
result_df <- data.frame(Gene = selected_genes, Occurrence = as.integer(gene_freq[selected_genes]))
gene_lists <- list()
for (fold in folds) {
  markers_rec <- read.csv(paste0("../survival/makers_rec_rank_fold", fold, ".csv"))
  gene_lists[[fold]] <- markers_rec$x
}

# Find genes present in all folds (intersection)
common_genes <- Reduce(intersect, gene_lists)
result_df <- result_df[result_df$Gene %in% common_genes,]
# Save result table
write.csv(result_df, "../survival/genes_occurrence_at_least_3.csv", row.names = FALSE)
