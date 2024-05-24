# 2019-11-21
# wgcna analysis to create signed correlation network

library(tidyverse)
library(WGCNA)

pheno <- read.csv("results/deseq2_moc_annot_table.csv", row.names = 1)
rownames(pheno) <- make.names(rownames(pheno))
exprs <- read.csv("results/deseq2_moc_counts_els_norm.csv", row.names = 1)


# Split into high and low risk subsets ------------------------------------


# labels
nSets <- 2
setLabels <- c("MOC", "BEN")

multiExpr <- vector(mode = "list", length = nSets)

multiExpr[[1]] <- list(data = t(exprs[,  rownames(pheno[pheno$Classification == "MOC", ])]))
multiExpr[[2]] <- list(data = t(exprs[,  rownames(pheno[pheno$Classification == "BEN", ])]))

exprSize <-  checkSets(multiExpr)

exprSize

gsg <- goodSamplesGenesMS(multiExpr, verbose = 3)
gsg$allOK


rm(exprs)

# cluster samples on euclid distance

sampleTrees <- list()
for(set in 1:nSets){
  sampleTrees[[set]] <- hclust(dist(multiExpr[[set]]$data), method = "average")
}

# pdf("results/plots/moc-wgcna-sample-clustering.pdf", width = 12, height = 12)
par(mfrow = c(2,1), mar = c(0, 4, 2, 0))
for(set in 1:nSets){
  plot(
    sampleTrees[[set]],
    main = paste("sample clustering on all genes in", setLabels[set]),
    xlab = "",cex = 0.7
  )
}


# Select soft-thresholding power ------------------------------------------

enableWGCNAThreads(nThreads = 4)

powers <- c(4:10, seq(12, 20, by = 2))

powerTables <- vector(mode = "list", length(nSets))

for(set in 1:nSets){
  powerTables[[set]] <- list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector = powers, verbose = 2)[[2]])
}

# plot results
colours <- c("black", "red")
plotCols <- c(2, 5, 6, 7)
colNames <- c(
  "Scale Free Topology Model Fit",
  "Mean connectivity",
  "Median connectivity",
  "Max connectivity"
)

yLim <- matrix(NA, nrow = 2, ncol = 4)
for(set in 1:nSets){
  for(col in 1:length(plotCols)){
    yLim[1, col] <- min(yLim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
    yLim[2, col] <- max(yLim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
  }
}

sizeGrWindow(8, 6)
par(mfrow = c(2,2))
par(mar = c(4.2, 4.2, 2.2, 0.5))
cex1 <- 0.7
# weird plot loop to look at topology fit results.
for(col in 1:length(plotCols)){
  for(set in 1:nSets){
    if(set == 1){
      plot(
        powerTables[[set]]$data[, 1],
        -sign(powerTables[[set]]$data[, 3]) * powerTables[[set]]$data[, 2],
        xlab = "Soft Threshold (power)", ylab = colNames[col], type = "n",
        ylim = yLim[, col], main = colNames[col]
      )
      addGrid()
    }
    if (col == 1){
      text(
        powerTables[[set]]$data[, 1],
        -sign(powerTables[[set]]$data[, 3]) * powerTables[[set]]$data[, 2],
        labels = powers, cex = cex1, col = colours[set]
      )
    } else {
      text(
        powerTables[[set]]$data[,1],
        powerTables[[set]]$data[,plotCols[col]],
        labels = powers,cex = cex1, col = colours[set]
      )
    }
    if(col == 1) {
      legend("bottomright", legend = setLabels, col = colours, pch = 20)
    } else{
      legend("topright", legend = setLabels, col = colours, pch = 20)
    }
  }
}

# from the above plot use a thresholding value of 10



# calculate network adjacency ---------------------------------------------

collectGarbage()

softPower <- 10
# nSets <- 1 testing
adjacencies <- array(0, dim = c(nSets, exprSize$nGenes, exprSize$nGenes))

for(set in 1:nSets){
  adjacencies[set, , ] <- abs(cor(multiExpr[[set]]$data, use = "p"))^softPower
}

collectGarbage()


hist(adjacencies[1, , ][upper.tri(adjacencies[1, , ])],ylim = c(0, 1E6))
hist(adjacencies[2, , ][upper.tri(adjacencies[1, , ])],ylim = c(0, 1E6))


# save adjacencies
genenames <- colnames(multiExpr[[1]]$data)
for(set in 1:nSets){
  write.table(
    adjacencies[set, , ],
    paste0("results/moc-wgcna-", softPower, "th-power-", setLabels[set], "-adjacency-matrix.txt"),
    col.names = genenames,
    row.names = genenames,
    quote = F,
    sep = "\t"
  )
}
