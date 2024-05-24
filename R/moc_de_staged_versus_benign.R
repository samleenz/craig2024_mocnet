# 2019-07-17
# Samuel Lee


# DE testing of benign versus invasive tumours
# using the MOC dataset

library(tidyverse)
library(DESeq2)
library(BiocParallel)
register(SnowParam(4))

args <- commandArgs(trailingOnly=TRUE)

datDir <- args[1]

annotRaw <- read_csv(file.path(datDir, "All survival_CN_Aug18.csv"))

table(annotRaw$Classification)

counts <- read_csv(file.path(datDir, "analysis_set_raw_counts_genenames.csv"))
counts <- counts[order(counts[[1]]), ]
countsSamples <- str_remove(colnames(counts)[-1], "GAMuT_")
countsSamples <- str_replace(countsSamples, "-", "/")


# double check there are no other name overlap issues
annot <- annotRaw %>%
  filter(
    GAMUT_ID %in% countsSamples,
    Classification %in% c("MOC", "BEN")
    ) %>%
  mutate(GAMUT_ID = paste0("x-", GAMUT_ID)) %>%
  mutate(combined = factor(paste(Classification, Stage))) %>%
  arrange(GAMUT_ID) %>%
  as.data.frame()

rownames(annot) <- annot$GAMUT_ID

glimpse(annot)

countsFilt <- as.data.frame(counts[, -1])
colnames(countsFilt) <- paste0("x-", countsSamples)
countsFilt <- countsFilt[! duplicated(counts[[1]]) ,] # drop duplicated gene rows
rownames(countsFilt) <- unique(counts[[1]])
countsFilt <- countsFilt[, annot$GAMUT_ID]




# Create DESeq2 object

dds <- DESeqDataSetFromMatrix(countsFilt, annot, ~ Classification)
# make comparison levels explicit
dds$Classification <- relevel(dds$Classification, ref = "BEN")


# run DESeq
dds <- DESeq(dds, parallel = T)

res <- results(dds)


plotMA(res)



resLFC <- lfcShrink(dds, coef="Classification_MOC_vs_BEN", type="apeglm")

plotMA(resLFC)


plotCounts(dds, order(res$padj)[[2]], intgroup = "combined")


vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup = "Classification")


# save results

write_csv(results(dds, tidy = T), file.path("results", "deseq2_moc_staged_vs_benign.csv"))

write_csv(annot, file.path("results", "deseq2_moc_annot_table.csv"))

# save normalised counts

normCounts <- counts(dds, normalized = TRUE)
write.csv(normCounts, file.path("results", "deseq2_moc_counts_els_norm.csv"))
