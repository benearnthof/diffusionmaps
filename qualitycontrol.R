# Quality control
# Motivation: 
# Low quality libraries in scRNA-seq data can arise from a variety of sources 
# such as cell damage during dissociation or failure in library preparation
# (inefficient reverse transcription, PCR amplification etc.)
# => Cells with low total counts, few expressed genes and high mitochondrial or 
# spike-in proportions. 
# => Can lead to misleading results in downstream analyses such as:
# Formation of own distinct clusters
# They distort the characterization of population heterogenity during variance 
# estimation or PCA 
# They contain genes that appear to be strongly "upregulated" due to aggressive 
# scaling to normalize for by small library sizes. 
# Increased scaliing in low quality libraries transforms small counts for these 
# transcripts into large normalized expression values, resulting in apparent 
# upregulation compared to other cells. 

# To avoid, or at least mitigate, these problems we need to remove these cells 
# at the start of the analysis. This step is called QC (Quality control.)

# demonstration on a small scRNA-seq dataset from Lun et al 2017

library(scRNAseq)
sce.416b <- LunSpikeInData(which = "416b") 
sce.416b$block <- factor(sce.416b$block)

sce.416b

# for each cell we calculate QC metrics using the
# perCellQCMetrics() function from scater
# sum: total count for each cell
# detected: number of detected genes
# subsets_Mito_percent: percentage of reads mapped to mitochondrial transcripts
# altexps_ERCC_percent: percentage of reads mapped to ERCC transcripts

# Retrieving the mitochondrial transcripts using genomic locations included in
# the row-level annotation for the SingleCellExperiment.
location <- rowRanges(sce.416b)
is.mito <- any(seqnames(location)=="MT")

# ALTERNATIVELY: using resources in AnnotationHub to retrieve chromosomal
# locations given the Ensembl IDs. Same result.
library(AnnotationHub)
ens.mm.v97 <- AnnotationHub()[["AH73905"]]
chr.loc <- mapIds(ens.mm.v97, keys=rownames(sce.416b),
                  keytype="GENEID", column="SEQNAME")
is.mito <- which(chr.loc=="MT")

library(scater)
df <- perCellQCMetrics(sce.416b, subsets=list(Mito=is.mito))
df

sce.416b <- addPerCellQC(sce.416b, subsets=list(Mito=is.mito))
colnames(colData(sce.416b))

# key assumption: QC metrics are independent of the biological state of each cell
# Poor values are presumed to be driven by technical factors rather than biological 
# processes, meaning that the subsequent removal of cells will not misrepresent the 
# biology in downstream analyses. Major violations of this assumption would 
# potentially result in the loss of cell types that have, say systematically low RNA
# content or high numbers of mitochondria. 

## Identifying low-quality cells

# library sizes below 100k reads, epxress fewer than 5000 genes
# have spike in proportions above 10 percent, or have mito proportions > 10%

qc.lib <- df$sum < 1e5
qc.nexprs <- df$detected < 5e3
qc.spike <- df$altexps_ERCC_percent > 10
qc.mito <- df$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito

# Summarize the number of cells removed for each reason.
DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
          SpikeProp=sum(qc.spike), MitoProp=sum(qc.mito), Total=sum(discard))

# these methods must be applied with care, depending on the type of sequencing 
# performed one must vary the quality thresholds appropriately

# identifying outliers with adaptive thresholds
# based on 3 median absolute deviations
# for this dataset we identify cells with log transformed library sizes
# that are more than 3MADs below the median. 
# Logtransform is used here to improve resolution at small values when type = "lower"
# It also guarantees that the threshold is a positive value
# Also it skews heavy right tailed distributions to be more "normal" looking

qc.lib2 <- isOutlier(df$sum, log=TRUE, type="lower")
# 4 outliers

# outlier identification for proportion based metrics are performed in a similar
# manner. Here no transformation is performed as we are aiming to identify and 
# remove large outliers that should already be clearly distinguishable from zero. 

qc.spike2 <- isOutlier(df$altexps_ERCC_percent, type="higher")
attr(qc.spike2, "thresholds")

qc.mito2 <- isOutlier(df$subsets_Mito_percent, type="higher")
attr(qc.mito2, "thresholds")

qc.nexprs2 <- isOutlier(df$detected, log=TRUE, type="lower")

discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

sum(discard2)
# only 6 outliers in total, this method is far more conservative

discard2 <- qc.lib2 | qc.nexprs2 | qc.spike2 | qc.mito2

# Summarize the number of cells removed for each reason.
DataFrame(LibSize=sum(qc.lib2), NExprs=sum(qc.nexprs2),
          SpikeProp=sum(qc.spike2), MitoProp=sum(qc.mito2), Total=sum(discard2))

# alternatively this can be performed in a single operation like so
reasons <- quickPerCellQC(df, percent_subsets=c("subsets_Mito_percent",
                                                "altexps_ERCC_percent"))
colSums(as.matrix(reasons))

# make sure to do these operations separately for every batch 
# but that assumes that the batches are all of similar quality

library(scRNAseq)
sce.grun <- GrunPancreasData()
sce.grun <- addPerCellQC(sce.grun)

# First attempt with batch-specific thresholds.
discard.ercc <- isOutlier(sce.grun$altexps_ERCC_percent,
                          type="higher", batch=sce.grun$donor)
with.blocking <- plotColData(sce.grun, x="donor", y="altexps_ERCC_percent",
                             colour_by=I(discard.ercc))

# Second attempt, sharing information across batches
# to avoid dramatically different thresholds for unusual batches.
discard.ercc2 <- isOutlier(sce.grun$altexps_ERCC_percent,
                           type="higher", batch=sce.grun$donor,
                           subset=sce.grun$donor %in% c("D17", "D2", "D7"))
without.blocking <- plotColData(sce.grun, x="donor", y="altexps_ERCC_percent",
                                colour_by=I(discard.ercc2))

gridExtra::grid.arrange(with.blocking, without.blocking, ncol=2)

# look for outlier batches first
ercc.thresholds <- attr(discard.ercc, "thresholds")["higher",]
ercc.thresholds
names(ercc.thresholds)[isOutlier(ercc.thresholds, type="higher")]

# check diagnostic plots

colData(sce.416b) <- cbind(colData(sce.416b), df)
sce.416b$block <- factor(sce.416b$block)
sce.416b$phenotype <- ifelse(grepl("induced", sce.416b$phenotype),
                             "induced", "wild type")
sce.416b$discard <- reasons$discard

gridExtra::grid.arrange(
  plotColData(sce.416b, x="block", y="sum", colour_by="discard",
              other_fields="phenotype") + facet_wrap(~phenotype) + 
    scale_y_log10() + ggtitle("Total count"),
  plotColData(sce.416b, x="block", y="detected", colour_by="discard", 
              other_fields="phenotype") + facet_wrap(~phenotype) + 
    scale_y_log10() + ggtitle("Detected features"),
  plotColData(sce.416b, x="block", y="subsets_Mito_percent", 
              colour_by="discard", other_fields="phenotype") + 
    facet_wrap(~phenotype) + ggtitle("Mito percent"),
  plotColData(sce.416b, x="block", y="altexps_ERCC_percent", 
              colour_by="discard", other_fields="phenotype") + 
    facet_wrap(~phenotype) + ggtitle("ERCC percent"),
  ncol=1
)

# Another useful diagnostic involves plotting the proportion of mitochondrial 
# counts against some of the other QC metrics. The aim is to confirm that there
# are no cells with both large total counts and large mitochondrial counts, 
# to ensure that we are not inadvertently removing high-quality cells that happen
# to be highly metabolically active (e.g., hepatocytes). In this case, we do not
# observe any points in the top-right corner in Figure 6.3.

plotColData(sce.416b, x="sum", y="subsets_Mito_percent", 
            colour_by="discard", other_fields=c("block", "phenotype")) +
  facet_grid(block~phenotype) +
  theme(panel.border = element_rect(color = "grey"))

# Comparison of the ERCC and mitochondrial percentages can also be 
# informative (Figure 6.4). Low-quality cells with small mitochondrial
# percentages, large spike-in percentages and small library sizes are 
# likely to be stripped nuclei, i.e., they have been so extensively 
# damaged that they have lost all cytoplasmic content. Conversely, cells
# with high mitochondrial percentages and low ERCC percentages may 
# represent undamaged cells that are metabolically active.

plotColData(sce.416b, x="altexps_ERCC_percent", y="subsets_Mito_percent",
            colour_by="discard", other_fields=c("block", "phenotype")) +
  facet_grid(block~phenotype) + 
  theme(panel.border = element_rect(color = "grey"))

# 6.5 Cell calling for Droplet Data
# droplet data may contain almost empty droplets with leftover RNA
# these need to be filtered out before the analysis.
library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
raw.path <- bfcrpath(bfc, file.path("http://cf.10xgenomics.com/samples",
                                    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"))
untar(raw.path, exdir=file.path(tempdir(), "pbmc4k"))
BiocManager::install("DropletUtils")
library(DropletUtils)
library(Matrix)
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names=TRUE)
sce.pbmc

# sharp transition between barcodes with large and small library sizes
bcrank <- barcodeRanks(counts(sce.pbmc))

# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)

abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Total UMI (Unique Molecular Identifier) count for each barcode in the 
# PBMC dataset, plotted against its rank (in decreasing order of total 
# counts). The inferred locations of the inflection and knee points 
# are also shown. 

# test for significant differences from empty droplets 
# emptyDrops performs Monte Carlo simulations to compute p-values,
# so we need to set the seed to obtain reproducible results.
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))

# See ?emptyDrops for an explanation of why there are NA values.
summary(e.out$FDR <= 0.001)

table(Sig=e.out$FDR <= 0.001, Limited=e.out$Limited)

# distribution of p values for droplet test should be close to uniform
# if there are spikes near zero we need to decrease "Lower"

set.seed(100)
limit <- 100   
all.out <- emptyDrops(counts(sce.pbmc), lower=limit, test.ambient=TRUE)
hist(all.out$PValue[all.out$Total <= limit & all.out$Total > 0],
     xlab="P-value", main="", col="grey80") 

# if we are satisfied with the results we can subset the data accordingly
sce.pbmc <- sce.pbmc[,which(e.out$FDR <= 0.001)]
