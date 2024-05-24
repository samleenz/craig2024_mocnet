# 20200203
# Samuel Lee

# Differential node analysis of WGCNA coexpression graphs
# 

library(nedd)
library(igraph)
library(tidyverse)
library(PCSF)

loadGraphToUniprot <- function(path){
  
  df <- readr::read_tsv(path) %>%
    dplyr::mutate(node_2 = stringr::str_remove(node_2, "^X")) 
  
  uniprot <- c(df$node_1, df$node_2) %>%
    unique() %>%
    gprofiler2::gconvert(target = "UNIPROTSWISSPROT") %>%
    dplyr::arrange(input_number) %>%
    dplyr::distinct(input_number, .keep_all = TRUE) %>%
    dplyr::select(input, target)
  
  df_uniprot <- df %>%
    dplyr::right_join(uniprot, by = c("node_1" = "input")) %>%
    dplyr::select(- node_1) %>%
    dplyr::rename(node_1 = target ) %>%
    dplyr::right_join(uniprot, by = c("node_2" = "input")) %>%
    dplyr::select(- node_2) %>%
    dplyr::rename(node_2 = target ) %>%
    tidyr::drop_na() %>%
    dplyr::select(node_1, node_2, weight)
  
  g_uni <- df_uniprot %>%
    arrange(node_1, node_2) %>% # this is a way to try make the vertex sequences match
    igraph::graph_from_data_frame(directed = F)# %>%
    # igraph::intersection(ascher_graph)
  
  return(g_uni)
}


# load wgcna weights, and make condition specific graphs
g_moc <-loadGraphToUniprot("results/moc-wgcna-10th-power-MOC-adjacency-3column-filterDot01.txt")

g_ben <- loadGraphToUniprot("results/moc-wgcna-10th-power-BEN-adjacency-3column-filterDot01.txt")

# drop nodes and edges that are not in the intersection of the graphs
# probably only need this once (could use E(g)$weight_2) but this is easier to make sense of.. 
g_moc_int <- intersection(g_moc, g_ben, keep.all.vertices = FALSE)
g_ben_int <- intersection(g_ben, g_moc, keep.all.vertices = FALSE)

# get differential edge scores then intersect with the ascher graph
g_diff <- nedd::diff_i(g_moc_int, g_ben_int, name1 = "weight_1", name2 = "weight_1", nameOut = "weight_diff")
g_diff_asch <- intersection(g_diff, ascher_graph) %>%
  delete_vertices(., V(.)[degree(.) == 0])

vcount(g_diff_asch); ecount(g_diff_asch)


# generate costs from differential weights
# cost is called $weight because of PCSF requirements
hist(E(g_diff_asch)$weight_diff)
hist(1 - abs(E(g_diff_asch)$weight_diff)) # making assumption that +ve/-ve diff_weights are of equal info
E(g_diff_asch)$weight <- (1 - abs(E(g_diff_asch)$weight_diff))
E(g_diff_asch)$weight <- ifelse(E(g_diff_asch)$weight < 0.1, 0.1, E(g_diff_asch)$weight)
hist(E(g_diff_asch)$weight)


###
# load SNP scores 
snp_prize <- read_tsv("results/snp_min_consequence_moc_noBen_scores.tsv", col_names = FALSE) %>%
  dplyr::filter(X2 > 0) %>%
  dplyr::filter(X1 %in% igraph::V(g_diff_asch)$name) %>%
  # dplyr::mutate(X2 = ifelse(X2 > 3, 0, X2)) %>% %>%
  dplyr::mutate(X2 = 4 - as.numeric(X2)) %>% # this is to flip the scores as 1 = most important
  tibble::deframe()

pcsf_res <- PCSF_rand(g_diff_asch, snp_prize, w = 2, b = 1, mu = 0.0005, n = 50)

vcount(pcsf_res)

pcsf_res_stats <- nedd::netStats(pcsf_res, weights = NA) %>%
  dplyr::mutate(protein = igraph::V(pcsf_res)$name) %>%
  dplyr::arrange(desc(betweenness), desc(drug_score))

pcsf_ranks <- pcsf_res_stats %>%
  dplyr::select(protein, drug_score, betweenness) %>%
  nedd::rankTable(nameCol = "protein", aggFUN = prod) %>%
  arrange(desc(aggRank))


save(g_diff_asch, pcsf_res, pcsf_res_stats, pcsf_ranks, file = "results/2020-02-06_moc_pcsf_wgcna_snp_results_filterdot01_diffI.RData")


