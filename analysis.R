options(stringsAsFactors = FALSE)
load("~/Desktop/Bioinformatics/Richard + Arif's project/Fluva_DP_PSet.RData")
fluvastatin.ic50 <- summarizeSensitivityProfiles(Fluva_DP_PSet, "IC50", drugs = "fluvastatin")
rnaseq <- summarizeMolecularProfiles(Fluva_DP_PSet, "rnaseq.counts", fill.missing = FALSE)
pr_cd_rows <- Biobase::fData(rnaseq)[,"gene_type"] == "protein_coding"

rnaseq <- rnaseq[pr_cd_rows,]

cell_lines <- rownames(drug.combo.data[,"5",])
col_val <- colnames(drug.combo.data[,"5",])

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

row_names_data <- c()
##Create vector which doesn't include duplicates
for (each in cell_lines){
  num<-nchar(each)
  fin <- substr(each, num-1, num-1)
  if (!fin=="_")row_names_data <- c(row_names_data, each)
}

##create updated matrix
updated_data = matrix(nrow = 44, ncol = 10, dimnames = list(row_names_data, col_val))
##find average of values and plug into new matrix
for (name in 1:length(row_names_data)){
  tot <- 0
  freq <- 0
  for (col in 1:ncol(drug.combo.data[,"5",])){
    for (row in 1:nrow(drug.combo.data[,"5",])){
      row_name <- cell_lines[row]
      if (!is.na(pmatch(row_names_data[name],row_name))){
        tot <- tot + drug.combo.data[,"5",][row,col]
        freq <- freq + 1
      }
    }
    val <- tot/freq
    updated_data[name, col] <- val
  }
}
sort(cellNames(Fluva_DP_PSet))
read.csv(file = "~/Desktop/cell_lines.csv")
d <- read.csv(file = "~/Desktop/cell_lines.csv")
row_names_data <- d[,3]
rownames(updated_data) <- row_names_data

##compute AUC Values
AUC_values <- c()

for (idx in 1:length(row_names_data)){
  auc_val <- computeAUC(concentration = names(updated_data[idx,]), viability = updated_data[idx,], conc_as_log = FALSE, viability_as_pct = FALSE)
  AUC_values <- c(AUC_values, auc_val)
}
AUC_values_fluva <- summarizeSensitivityProfiles(Fluva_DP_PSet, "AUC", fill.missing = FALSE)[2,]
names(AUC_values) <- row_names_data


## Waterfall to determine resistant vs. sensitive
sens_res_values <- PharmacoGx:::callingWaterfall(AUC_values, type = "AUC", intermediate.fold = 1)
sens_res_values_fluva <- PharmacoGx:::callingWaterfall(AUC_values_fluva, type = "AUC", intermediate.fold = 1)
sens_res_values_fluva <- sort(sens_res_values_fluva)
source("~/Desktop/callingWaterfall.R")
plot_sens_res_values_fluva <- callingWaterfall(AUC_values_fluva, type = "AUC", intermediate.fold = 1, plot = TRUE)
##update rna column values
##deseq
library(Biobase)
exprs(rnaseq) <- round(2^exprs(rnaseq) - 1)
library(SummarizedExperiment)
library(DESeq2)

RNASeqSE <- makeSummarizedExperimentFromExpressionSet(rnaseq)
RNASeqSE <- RNASeqSE[,names(sens_res_values)]
##rerun from here
colData(RNASeqSE)$fluva_dp <- sens_res_values[colnames(RNASeqSE)]
# colData(RNASeqSE)$fluva <- sens_res_values_fluva[colnames(RNASeqSE)]
dds <- DESeqDataSet(RNASeqSE, design = ~ fluva_dp)
dds <- DESeq(dds)

##results
res <- results(dds, alpha = 0.05, cooksCutoff= 4/42)
res$Symbol <- rowData(RNASeqSE)$Symbol[match(rownames(res) ,rowData(RNASeqSE)$gene_id)]
resSig <- subset(res, res$padj < 0.05)
resSig <- resSig[order(resSig$log2FoldChange),]
# rownames(resSig) <- NULL
write.csv(resSig,"resSig.csv")


RNASeqSE <- makeSummarizedExperimentFromExpressionSet(rnaseq)
RNASeqSE <- RNASeqSE[,intersect(colnames(RNASeqSE),names(sens_res_values_fluva))]
##rerun from here
colData(RNASeqSE)$fluva <- sens_res_values_fluva[colnames(RNASeqSE)]
##deseq for fluva
dds_fluva <- DESeqDataSet(RNASeqSE, design = ~ fluva)
dds_fluva <- DESeq(dds_fluva)

##results for fluva
res <- results(dds_fluva, alpha = 0.05, cooksCutoff= 4/42)
res$Symbol <- rowData(RNASeqSE)$Symbol[match(rownames(res) ,rowData(RNASeqSE)$gene_id)]
resSig <- subset(res, res$padj < 0.05)
resSig <- resSig[order(resSig$log2FoldChange),]
write.csv(resSig,"resSig_fluva2.csv")

dose_response_data <- Fluva_DP_PSet@sensitivity$raw[c("CAMA-1_fluvastatin_rep_1","MDA-MB-231_fluvastatin_rep_1"),,]
concentrations <- dose_response_data[,,"Dose"]
viabilities <- dose_response_data[,,"Viability"]

concentrations <- list("CAMA-1" = concentrations[1,], "MDA-MB-231" =concentrations[2,])
viabilities <- list("CAMA-1" = viabilities[1,], "MDA-MB-231" =viabilities[2,])

drugDoseResponseCurve(concentrations = concentrations, viabilities = viabilities, legends.label = "auc_recomputed", plot.type = "Fitted")

RNASeqSEsig <- RNASeqSE[rownames(resSig),]
assay(RNASeqSEsig)

log2rna <- log2(assay(RNASeqSEsig)+1)
log2rna_z <- t(apply(log2rna, 1, scale))
colnames(log2rna_z) <- colnames(RNASeqSEsig)

library(pheatmap)
pheatmap(log2rna_z[,names(sens_res_values_fluva)], labels_row = NA,annotation_col = data.frame(sens_res_values_fluva), cluster_rows = FALSE, cluster_cols = FALSE)
