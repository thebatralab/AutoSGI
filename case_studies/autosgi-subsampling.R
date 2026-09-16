# this script checks for reproducibility of neuropathological significance under subsampling.
set.seed(7890)
# ---- Libraries --------------------------------------------------------------
library(cluster)   # for daisy/dist if needed
library(tidyverse)
library(magrittr)
library(maplet)
library(autosgi)
library(openxlsx)


# ---- Custom functions (provides ordinal_test) -------------------------------
source("/case_studies/0_custom_functions.R")

# ---- Input data (ROSMAP, scaled) --------------------------------------------
load('/case_studies/rosmap-results/rosmap-metabo-autosgi-data.rds')
load('/case_studies/rosmap-results/rosmap-metabo-hierarchical-selection-358_sgi_data.rds')

# ---- Seven-metabolite signature ---------------------------------------------
seven_metabolites <- c(
  "glutamate", "N-acetylglycine", "2-aminoadipate", "guanidinoacetate",
  "glycerophosphoethanolamine", "X - 25020", "glycerophosphorylcholine (GPC)"
)

# feature selection (subset to the signature)
metabolite_matrix <- metabolite_matrix [, which(rowData(rosmap_data)$name %in% seven_metabolites)]
d <- metabolite_matrix %>% as.matrix()
#colnames(d) <- rownames(phenotype_data)
sample_composition


## ---- INPUTS YOU PROVIDE ----
#matrix/data.frame, rows = individuals, cols = the 7 fixed signature metabolites
X <- d
# named vector of original cluster labels for each individual
#                (names = individual IDs, matching rownames(X))
orig_clusters <- as.vector(sample_composition$l2)-1
names(orig_clusters) <- sample_composition$sampleids
#number of clusters to cut the tree into (match your original)
k <- 2
#number of subsamples (e.g., 1000)
n_iter <- 1000
subsample_frac <- 0.8
results <- vector("list", n_iter <-1000)
cogdx_results <- vector("list", n_iter <-1000)
braaksc_results <- vector("list", n_iter <-1000)
ceradsc_results <- vector("list", n_iter <-1000)
subsample_frac <- 0.8

## ---- CLUSTERING FUNCTION (adapt to your pipeline) ----
# Replace this with however AutoSGI/SGI clusters samples on the signature.
# Here: Euclidean + Ward + cut into k groups, to match your paper.
cluster_fn <- function(Xsub, k) {
  d  <- dist(Xsub, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  cutree(hc, k = k)
}

## ---- JACCARD BETWEEN TWO CLUSTERINGS (per original cluster) ----
# For each original cluster, find the best-matching subsample cluster
# and record the Jaccard overlap (computed only over the shared individuals).
jaccard_per_cluster <- function(orig, sub, ids) {
  # orig, sub: named cluster-label vectors; ids: individuals in the subsample
  orig <- orig[ids]
  sub  <- sub[ids]
  orig_labels <- unique(orig)
  sapply(orig_labels, function(oc) {
    orig_members <- names(orig)[orig == oc]
    # best-matching subsample cluster
    best <- max(sapply(unique(sub), function(sc) {
      sub_members <- names(sub)[sub == sc]
      inter <- length(intersect(orig_members, sub_members))
      union <- length(union(orig_members, sub_members))
      if (union == 0) 0 else res <- inter / union
      return(res)
    }))
    best
  })
}

## ---- MAIN LOOP ----
set.seed(1)
n <- nrow(X)
ids_all <- rownames(X)
k <- length(unique(orig_clusters))

for (i in seq_len(n_iter)) {
  # subsample WITHOUT replacement, fixed fraction
  idx <- sample(ids_all, size = floor(subsample_frac * n), replace = FALSE)
  Xsub <- X[idx, , drop = FALSE]
  Xpheno <- phenotype_data[idx, , drop = FALSE]
  sub_clusters <- cluster_fn(Xsub, k)
  names(sub_clusters) <- idx
  # Jaccard of each original cluster vs best subsample match, over shared ids
  results[[i]] <- jaccard_per_cluster(orig_clusters, sub_clusters, ids = idx)
  cogdx_results[[i]] <- ordinal_test (Xpheno$cogdx, sub_clusters) %>% unlist()
  braaksc_results[[i]] <- ordinal_test (Xpheno$braaksc, sub_clusters) %>% unlist()
  ceradsc_results[[i]] <- ordinal_test (Xpheno$ceradsc, sub_clusters) %>% unlist()
}

# results formatting

jaccard_df <- do.call(rbind, results) %>% data.frame() 
cogdx_df <- do.call(rbind, cogdx_results) %>% data.frame() 
cogdx_df %<>% mutate(p_adj=p.adjust(pval, method="bonferroni")) %>% data.frame()
braaksc_df <- do.call(rbind, braaksc_results) %>% data.frame()
braaksc_df %<>%mutate(p_adj=p.adjust(pval, method="bonferroni")) %>% data.frame()
ceradsc_df <- do.call(rbind, ceradsc_results) %>% data.frame()
ceradsc_df %<>% mutate(p_adj=p.adjust(pval, method="bonferroni")) %>% data.frame()

# writing out results
outcomes <- c("jaccard", "cogdx", "braaksc", "ceradsc")
wb <- openxlsx::createWorkbook()
# create worksheet
openxlsx::addWorksheet(wb,"jaccard_index")
# write data
openxlsx::writeData(wb, "jaccard_index", jaccard_df,rowNames = F, colNames = T)
# create worksheet
openxlsx::addWorksheet(wb,"cogdx_stats")
# write data
openxlsx::writeData(wb, "cogdx_stats", cogdx_df,rowNames = F, colNames = T)
# create worksheet
openxlsx::addWorksheet(wb,"braaksc_stats")
# write data
openxlsx::writeData(wb, "braaksc_stats", braaksc_df,rowNames = F, colNames = T)
# create worksheet
openxlsx::addWorksheet(wb,"ceradsc_stats")
# write data
openxlsx::writeData(wb, "ceradsc_stats", ceradsc_df,rowNames = F, colNames = T)
# write workbook
openxlsx::saveWorkbook (wb, file="/case_studies/autosgi_subsampling.xlsx", overwrite=TRUE)
