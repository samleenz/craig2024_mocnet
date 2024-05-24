

# load packages
library(tidyverse)
library(nedd)
library(igraph)
library(ggraph)
library(patchwork)

# load the network analysis results 
load("data/2020-02-06_moc_pcsf_wgcna_snp_results_filterdot01_diffI.RData")


# Calculate betweenness centrality
g_diff_btwn <- betweenness(g_diff_asch)

# Generate network statistics and combine with betweenness centrality
networkStats <- pcsf_res |>
  netStats() |>
  rownames_to_column(var = "uniprot") |>
  mutate(betweenness = g_diff_btwn[uniprot]) |>
  dplyr::select(uniprot, betweenness, drug_score) |>
  mutate(occurences = V(pcsf_res)$prize)

# Create a results table with rank product
results_table <- left_join(
  networkStats,
  rankTable(networkStats, nameCol = "uniprot", aggFUN = prod) |>
    dplyr::select(uniprot, aggRank)
) |>
  rename(rank_product = aggRank) |>
  arrange(desc(rank_product))

# Extract the top 10 Uniprot IDs based on rank product
uniprot_top10 <- results_table |>
  slice_max(rank_product, n = 10) |>
  pull(uniprot)

# Calculate and display the correlation between drug score and betweenness
networkStats |> 
  summarise(cor = cor(drug_score, betweenness, method = "pearson"))


cor.test(networkStats$drug_score, networkStats$betweenness, method = "pearson")


# Plot drug score versus betweenness and
networkStats |> 
  ggplot(aes(drug_score, log2(betweenness))) +
  geom_point() +
  ggpubr::theme_pubr(border = TRUE) 


## Make plot of largest 5 clusters
pcsf_clusters <- clusters(pcsf_res)
V(pcsf_res)$cluster <- pcsf_clusters$membership

pcsf_clusters$csize |> enframe("Cluster", "Size")

big5 <- pcsf_clusters$csize |> enframe("Cluster", "Size") |>
  # slice_max(Size, n = 5)
  filter(Size > 2)

pcsf_clusters$membership[pcsf_clusters$membership %in% big5$Cluster]

# Prepare the network to plot

graph_to_plot <- nedd::getSubnet(pcsf_res, (V(pcsf_res)[cluster %in% big5$Cluster])$name)
graph_to_plot <- delete_edge_attr(graph_to_plot, "weight")

V(graph_to_plot)$drug_score <- networkStats |>
  mutate(cluster = pcsf_clusters$membership) |>
  filter(cluster %in% big5$Cluster) |> 
  pull(drug_score)

V(graph_to_plot)$betwixt <- networkStats |>
  mutate(cluster = pcsf_clusters$membership) |>
  filter(cluster %in% big5$Cluster) |>
  pull(betweenness)

  
V(graph_to_plot)$evidence <- ifelse(V(graph_to_plot)$type == "Steiner", "Indirect", "Direct")

gene_names <- gprofiler2::gconvert(results_table$uniprot, filter_na = FALSE) |>
  group_by(input) |>
  slice_head(n = 1) |>
  ungroup() |>
  select(input, name) |>
  deframe()

V(graph_to_plot)$label <- gene_names[V(graph_to_plot)$name]



# Generate network diagrams -----------------------------------------------

  
  symbolToPlot <- c("CCNA2", "TRIP13", "CDK1", "CDC20", "PRC1")
  uniprotToPlot <- c("P20248", "Q15645", "P06493", "Q12834", "O43663")
  
  ##
  p1 <- ggraph(graph_to_plot, layout = "kk") +
    geom_edge_link() +
    geom_node_point(
      aes(size = drug_score, fill = log2(betwixt + 1), shape = evidence),
      colour = "white",
      stroke = 2
      ) +
    geom_node_text(aes(label = label, filter = label %in% symbolToPlot), repel = TRUE) +
    theme_void() +
    # scale_fill_viridis_c() +
    scico::scale_fill_scico(palette = "batlow", direction = 1) +
    scale_shape_manual(values = c(21,23)) +
    guides(
      shape = guide_legend(override.aes = list(color = "black", fill = "black", size = 3)),
      size = guide_legend(override.aes = list(color = "black", fill = "black"))
      )+
    labs(
      shape = "Evidence",
      fill = "log2(BC + 1)",
      size = "Drug score"
    )
  
  p2 <- ggraph(graph_to_plot, layout = "kk") +
    geom_edge_link() +
    geom_node_point(
      aes(size = drug_score, fill = log2(betwixt + 1), shape = evidence),
      colour = "white",
      stroke = 2
    ) +
    # geom_node_text(aes(label = label, filter = label %in% symbolToPlot), repel = TRUE) +
    theme_void() +
    # scale_fill_viridis_c() +
    scico::scale_fill_scico(palette = "batlow", direction = 1) +
    scale_shape_manual(values = c(21,23)) +
    guides(
      shape = guide_legend(override.aes = list(color = "black", fill = "black", size = 3)),
      size = guide_legend(override.aes = list(color = "black", fill = "black"))
    )+
    labs(
      shape = "Evidence",
      fill = "log2(BC + 1)",
      size = "Drug score"
    ); p2
  
  

# Save the plots and the results table -------------------------------------


  
  ggsave(
    filename = "plots/pcsf_subnetwork_figure_bc_drug_evidence.png",
    plot = p1
    )
  
  ggsave(
    filename = "plots/pcsf_subnetwork_figure_bc_drug_evidence_noText.png",
    plot = p2
  )
  

  left_join(
    results_table,
    igraph::as_data_frame(pcsf_res, what = "vertices")[,c(1, 3, 4)],
    by = c("uniprot" = "name")
  ) |>
    mutate(across(c("drug_score", "rank_product"), ~ round(., digits = 2))) |>
    mutate(Evidence = ifelse(type == "Steiner", "Indirect", "Direct")) |>
    mutate(Symbol = gene_names[uniprot]) |>
    select(uniprot, Symbol, everything(), -type) |>
    rename(
      `Uniprot ID` = uniprot,
      `HGNC Symbol` = Symbol,
      Betweenness = betweenness,
      `Drug score` = drug_score,
      Occurences = occurences,
      `Rank product` = rank_product
    ) |>
    write_csv("results/PCSF_rankproduct_table.csv")

  

  