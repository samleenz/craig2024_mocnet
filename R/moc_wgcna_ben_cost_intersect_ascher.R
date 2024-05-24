# 2020-01-09
# Samuel Lee

# intersect the wgcna cost graph with the ascher graph

library(tidyverse)
library(nedd)
library(igraph)

ben_cost <- readr::read_tsv(
  "results/wgcna-10th-power-BEN-edgelist-cost-filterDot01.txt"
)

ben_graph <- igraph::graph_from_data_frame(ben_cost, directed = F)

ben_ascher_graph <- igraph::intersection(ben_graph, ascher_graph, keep.all.vertices = FALSE)

ben_ascher_3col <- igraph::as_data_frame(ben_ascher_graph) %>%
  dplyr::rename(node_1 = from, node_2 = to)

readr::write_tsv(
  ben_ascher_3col,
  "results/wgcna-10th-power-BEN-edgelist-cost-filterDot01-ascherIntersect.txt"
)
