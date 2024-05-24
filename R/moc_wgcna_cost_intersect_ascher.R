# 2020-01-09
# Samuel Lee

# intersect the wgcna cost graph with the ascher graph

library(tidyverse)
library(nedd)
library(igraph)

moc_cost <- readr::read_tsv(
  "results/wgcna-10th-power-MOC-edgelist-cost-filterDot01.txt"
)

moc_graph <- igraph::graph_from_data_frame(moc_cost, directed = F)

moc_ascher_graph <- igraph::intersection(moc_graph, ascher_graph, keep.all.vertices = FALSE)

moc_ascher_3col <- igraph::as_data_frame(moc_ascher_graph) %>%
  dplyr::rename(node_1 = from, node_2 = to)

readr::write_tsv(
  moc_ascher_3col,
  "results/wgcna-10th-power-MOC-edgelist-cost-filterDot01-ascherIntersect.txt"
  )
