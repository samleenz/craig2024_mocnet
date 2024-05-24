# 2019-04-16
# Samuel Lee

#### General use functions
library(magrittr)



getFreqs <- function(vec){
  # return a named vector of frequencies for the 5 copy number states
  # for the MOC CNV data. input is a vector of CNV categories for a feature (gene)
  freqs <- c(
    "HCL" = sum(vec == "Homozygous Copy Loss"),
    "CNL" = sum(vec == "CN Loss"),
    "DIP" = sum(vec == "Diploid"),
    "CNG" = sum(vec == "CN Gain"),
    "HCG" = sum(vec == "High Copy Gain")
  )
  return(freqs)
}


softmax <- function(freqs, pc = 0){
  # returns the softmax log-prob of a vector of observed class frequencies
  # pc is a pseudo-count to add to all frequencies, no reason to use this
  sm <-  exp(freqs + pc) / sum(exp(freqs + pc))
  return(sm)
}

standardise <- function(x, ...){
  # scale a numeric vector to between [0,1]
  (x - min(x, ...)) / (max(x, ...) - min(x, ...))
}

getConsequence <- function(df, annot, cat, fun, minProp = 0){
  # extract relevant samples
  snpCAT <- df[df$Sample %in% annot$GAMUT_ID[annot$Classification == cat], ]
  # get popn size
  numCAT <- sum(annotRaw$Classification == cat, na.rm = TRUE)
  # get biggest snp effect per gene for each popn
  snpCATmax <- snpCAT %>%
    group_by(SYMBOL) %>%
    summarise(score = fun(Consequence_Rank), n = n()) %>%
    mutate(prop = n / numCAT) %>%
    arrange(desc(n)) %>%
    filter(prop > minProp)
  
  return(snpCATmax)
}

getScore <- function(snp, idx){
  # get uniprot names
  uniprot <- gprofiler2::gconvert(snp$SYMBOL, target = "UNIPROTSWISSPROT") %>%
    dplyr::arrange(input_number) %>%
    dplyr::distinct(input, .keep_all = TRUE)
  
  score <- dplyr::right_join(snp, dplyr::select(uniprot, input, target), by = c("SYMBOL" = "input")) %>%
    dplyr::select(target, score)
  
  score_full <- dplyr::left_join(idx, score, by = c("uniprotswissprot" = "target")) %>%
    tidyr::replace_na(list(score = 0)) %>%
    dplyr::select(-i)
  
  return(score_full)
}