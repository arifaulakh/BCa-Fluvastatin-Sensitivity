load("~/Desktop/Bioinformatics/CCLE.CTRPv2.RData")
features <- fNames(CCLE.CTRPv2, "rnaseq")
sig.rnaseq <- drugSensitivitySig(pSet = CCLE.CTRPv2, mDataType = "rnaseq", drugs = c("fluvastatin"),features=features)
sig.rnaseq_breast <- drugSensitivitySig(pSet = CCLE.CTRPv2, mDataType = "rnaseq", drugs = c("fluvastatin"),features=features,tissues = "breast")
sig.rnaseq <- sig.rnaseq[,1,]
sig.rnaseq_breast <- sig.rnaseq_breast[,1,]
symbols <- res[,7]
sig_symbols <- resSig[,7]
S <- symbols %in% sig_symbols
SIG <- res[S,]
ids <- rownames(SIG)
for (i in 1:length(ids)){
  temp <- strsplit(ids[i], split = ".",fixed = TRUE)[[1]]
  ids[i] <- temp[1]
}
sigrnaseq <- sig.rnaseq[intersect(ids, rownames(sig.rnaseq)),]
write.csv(sigrnaseq, "sigrnaseq.csv")
sigrnaseq_breast <- sig.rnaseq_breast[intersect(ids, rownames(sig.rnaseq_breast)),]
write.csv(sigrnaseq_breast, "sigrnaseq_breast.csv")
