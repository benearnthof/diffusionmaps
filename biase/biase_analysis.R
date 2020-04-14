# run create_sce.R first then set wd to this folder
### DATA
d <- read.table("data.txt", header = T)
d <- d[!duplicated(d$ID),]
rownames(d) <- d$ID
d <- d[,2:ncol(d)]

### ANNOTATIONS
ann <- read.table("biase_cell_types.txt", stringsAsFactors = F)

### SINGLECELLEXPERIMENT
source("../create_sce.R")
sceset <- create_sce_from_normcounts(d, ann)
# convert ensembl ids into gene names
# gene symbols will be stored in the feature_symbol column of fData
sceset <- getBMFeatureAnnos(
  sceset, filters="ensembl_gene_id",
  biomart="ensembl", dataset="mmusculus_gene_ensembl")
# remove features with duplicated names
sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]

saveRDS(sceset, "biase.rds")

biase <- readRDS("biase.rds")

library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
chr.loc <- mapIds(ens.mm.v97, keys=rownames(biase),
                  keytype="GENEID", column="SEQNAME")
is.mito <- which(chr.loc=="MT")

df <- perCellQCMetrics(biase, subsets=list(Mito=is.mito))
df

biase_qc <- addPerCellQC(biase, subsets=list(Mito=is.mito))
colnames(colData(biase_qc))

# identifying low quality gene reads
reasons <- quickPerCellQC(df, percent_subsets=c("subsets_Mito_percent"))
colSums(as.matrix(reasons))

# four cells will be discarded
# First attempt with batch-specific thresholds.
colnames(colData(biase_qc))
discard.mito <- isOutlier(biase_qc$subsets_Mito_percent, type="higher")
# only two would be discarded
plotColData(biase_qc, x="cell_type1", y="subsets_Mito_percent",
                             colour_by=I(discard.mito))

####### to understand needed structure for plotColData
# library(scRNAseq)
# sce.grun <- GrunPancreasData()
# sce.grun <- addPerCellQC(sce.grun)
# # First attempt with batch-specific thresholds.
# # is.outlier takes colnames of coldata
# # plotColData too
# discard.ercc <- isOutlier(sce.grun$altexps_ERCC_percent,
#                           type="higher", batch=sce.grun$donor)
# with.blocking <- plotColData(sce.grun, x="donor", y="altexps_ERCC_percent",
#                              colour_by=I(discard.ercc))

# checking diagnostic plots
biase_3 <- biase_qc

colData(biase_3) <- cbind(colData(biase_3), df)
colnames(colData(biase_3))
biase_3$Type <- factor(biase_3$cell_type1)
biase_3$Discard <- reasons$discard

plotColData(biase_3, x="Type", y="sum", colour_by="Discard") + 
  scale_y_log10() + 
  ggtitle("Total Count of Reads per Cell: Biase") +
  ylab("Total Count") +
  scale_fill_manual(values = c("#0000FF", "#ff0000")) +
  labs(fill = "Discard")

plotColData(biase_3, x="Type", y="detected", colour_by="Discard") + 
  scale_y_log10() + 
  ggtitle("Detected Features per Cell: Biase") +
  ylab("Detected Features") +
  scale_fill_manual(values = c("#0000FF", "#ff0000")) +
  labs(fill = "Discard")

plotColData(biase_3, x="Type", y="subsets_Mito_percent", 
            colour_by="Discard") +
  ggtitle("Percentage of Mitochondrial RNA: Biase") +
  ylab("Mitochondrial RNA in %") +
  scale_fill_manual(values = c("#0000FF", "#ff0000")) +
  labs(fill = "Discard")

# trying out droplet based diagnostics
library(DropletUtils)
library(Matrix)

bcrank <- barcodeRanks(counts(biase_3))
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

biase_3_filtered <- biase_3[,!reasons$discard]
biase_3 #56
biase_3_filtered #52
