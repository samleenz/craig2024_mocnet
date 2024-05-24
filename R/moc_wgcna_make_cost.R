# 2019-11-23
# Samuel Lee

# Convert WGCNA scores to edge costs 
# 

library(tidyverse)
library(nedd)
library(igraph)


moc_wgcna <- readr::read_tsv(
  "results/moc-wgcna-10th-power-MOC-adjacency-3column-filterDot01.txt"
)

# this is for the gene "7SK", which is made to be X7SK.. .
moc_wgcna <- moc_wgcna %>%
  mutate(node_2 = stringr::str_remove(node_2, "^X")) 


plot(density(moc_wgcna$weight))

# get unirpot IDs for gene names
uniprot <- c(moc_wgcna$node_1, moc_wgcna$node_2) %>%
  unique() %>%
  gprofiler2::gconvert(target = "UNIPROTSWISSPROT") %>%
  dplyr::arrange(input_number) %>%
  dplyr::distinct(input_number, .keep_all = TRUE) %>%
  dplyr::select(input, target)

# replace gene names with uniprot
# transform weights (corr) into cost ( 1 - weight) and 
# drop weights of 1 (self edges)
moc_wgcna_uni <- moc_wgcna %>%
  dplyr::right_join(uniprot, by = c("node_1" = "input")) %>%
  dplyr::select(- node_1) %>%
  dplyr::rename(node_1 = target ) %>%
  dplyr::right_join(uniprot, by = c("node_2" = "input")) %>%
  dplyr::select(- node_2) %>%
  dplyr::rename(node_2 = target ) %>%
  dplyr::mutate(cost = 1 - weight) %>%
  dplyr::filter(cost > 0) %>%
  dplyr::select(- weight) %>%
  dplyr::distinct(node_1, node_2, .keep_all = TRUE)
  

write_tsv(moc_wgcna_uni, "results/wgcna-10th-power-MOC-edgelist-cost-filterDot01.txt")

