# Samuel Lee
# make adjacency matrices from WGCNA into 3 column tables.
# Much more efficient to compute over


moc_mat <- as.matrix(read.table("results/moc-wgcna-10th-power-MOC-adjacency-matrix.txt", as.is = TRUE))
moc_long <- reshape2::melt(moc_mat)
colnames(moc_long) <- c("node_1", "node_2", "weight")
readr::write_tsv(moc_long, "results/moc-wgcna-10th-power-MOC-adjacency-3column.txt")

ben_mat <- as.matrix(read.table("results/moc-wgcna-10th-power-BEN-adjacency-matrix.txt", as.is = TRUE))
ben_long <- reshape2::melt(ben_mat)
colnames(ben_long) <- c("node_1", "node_2", "weight")
readr::write_tsv(ben_long, "results/moc-wgcna-10th-power-BEN-adjacency-3column.txt")
