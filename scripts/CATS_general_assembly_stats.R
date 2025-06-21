#!/usr/bin/env Rscript
#General transcriptome assembly statistics calculation script
#Loading the required package
suppressPackageStartupMessages(library(data.table))

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
THREAD_N <- as.numeric(ext_args[1])
TR_LEN_DT_PATH <- ext_args[2]
OUT_PREF <- ext_args[3]

#Importing transcript lengths from file
setDTthreads(THREAD_N)
tr_lengths <- fread(TR_LEN_DT_PATH, header = F, select = 2)
setnames(tr_lengths, "tr_len")
tr_lengths[, "tr_id" := .I]

#Calculating descriptive statistics of transcript length
tr_N <- tr_lengths[, .N]
base_N_sum <- sum(tr_lengths[, tr_len])
tr_len_summary <- summary(tr_lengths[, tr_len])

tr_longer_than_N <- c()
len_thr <- c(200, 500, 1000, 5000, 10000, 20000)
for (i in 1 : length(len_thr)) {
 tr_longer_than_N[i] <- tr_lengths[tr_len > len_thr[i], .N]
}

tr_longer_than_label <- c()
for (i in 1 : length(len_thr)) {
 tr_longer_than_label[i] <- paste("N, % transcripts longer than", len_thr[i], "bp")
}

#Calculating Nx and Lx metrics
setorder(tr_lengths, -tr_len)
tr_lengths[, "tr_len_cumsum" := cumsum(tr_len)]
N_L_50 <-  tr_lengths[tr_len_cumsum >= base_N_sum * 0.5][1]
N50 <- N_L_50[, tr_len]
L50 <- N_L_50[, tr_id]

N_L_90 <-  tr_lengths[tr_len_cumsum >= base_N_sum * 0.9][1]
N90 <- N_L_90[, tr_len]
L90 <- N_L_90[, tr_id]

#Writing general transcriptome assembly statistics to file
write.table(data.table(parameter = c("N transcripts", "Total assembly length (bp)", tr_longer_than_label, "Mean transcript length (bp)", "Median transcript length (bp)", "Transcript length IQR (bp)", "Transcript length range (bp)", "N50 (bp)", "L50", "N90 (bp)", "L90"), value = c(tr_N, base_N_sum, paste0(tr_longer_than_N, ", ", round(100 * tr_longer_than_N / tr_N, 2), "%"), round(tr_len_summary[4], 2), tr_len_summary[3], paste(tr_len_summary[2], tr_len_summary[5], sep = "-"), paste(tr_len_summary[1], tr_len_summary[6], sep = "-"), N50, L50, N90, L90)), file = paste(OUT_PREF, "general_statistics_table.tsv", sep = "_"), sep = "\t", row.names = F, col.names = F, quote = F)