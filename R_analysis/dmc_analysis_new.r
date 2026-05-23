library(pacman)
library(matrixTests)


p_load(data.table, bsseq, ggplot2, DSS, annotatr, GenomicRanges,
       TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db, pheatmap)

options(scipen=999)

outdatadir <- "/home/genetics/homework2/output/dmc_analysis/"
dir.create(outdatadir, recursive=TRUE)

sample_info <- data.frame(
  sample_id = c("SRR11647648", "SRR11647649", "SRR11647658", "SRR11647659"),
  condition = c("normal", "normal", "tumor", "tumor"),
  stringsAsFactors = FALSE
)

sample_colors <- c(
  "SRR11647648" = "#e89eb8",
  "SRR11647649" = "#332932",
  "SRR11647658" = "#bb9999",
  "SRR11647659" = "#bcbd22"
)


setDT(sample_info)


# ============================================================
# LOAD DATA AND CREATE BSSEQ OBJECT
# ============================================================
library(data.table)

files <- c(
  SRR11647648 = "~/homework2/methylation_chr20/SRR11647648_methylation_CpG.bedGraph",
  SRR11647649 = "~/homework2/methylation_chr20/SRR11647649_methylation_CpG.bedGraph",
  SRR11647658 = "~/homework2/methylation_chr20/SRR11647658_methylation_CpG.bedGraph",
  SRR11647659 = "~/homework2/methylation_chr20/SRR11647659_methylation_CpG.bedGraph"
)

sample_cols <- c("SRR11647648", "SRR11647649", "SRR11647658", "SRR11647659")


#--------------
ref <- fread(files[1], skip = 1, header = FALSE)

setnames(ref, c("chr","start","end","meth_pct","M","U"))

ref_id <- paste(ref$chr, ref$start, sep = "_")

M_list <- list()
Cov_list <- list()

for (s in names(files)) {

  dat <- fread(files[s], skip = 1, header = FALSE)
  setnames(dat, c("chr","start","end","meth_pct","M","U"))

  dat$id <- paste(dat$chr, dat$start, sep = "_")

  # align to reference CpGs
  dat <- dat[match(ref_id, dat$id)]

  # ensure numeric
  M <- as.numeric(dat$M)
  U <- as.numeric(dat$U)

  M[is.na(M)] <- 0
  U[is.na(U)] <- 0

  M_list[[s]] <- M
  Cov_list[[s]] <- M + U
}


M_matrix <- do.call(cbind, M_list)
Cov_matrix <- do.call(cbind, Cov_list)

colnames(M_matrix) <- names(files)
colnames(Cov_matrix) <- names(files)


stopifnot(all(M_matrix <= Cov_matrix, na.rm = TRUE))
#------------------


BS <- BSseq(
  chr = ref$chr,
  pos = ref$start,
  M = M_matrix,
  Cov = Cov_matrix,
  sampleNames = names(files)
)

cat("BSseq object created!\n")
cat("Dimensions:", dim(BS), "\n")


# ============================================================
# FILTER BY COVERAGE
# ============================================================
# Sukuriame methylation_dt ir coverage_dt iš BS objekto

methylation_dt <- data.table(
  cpg_id = paste(ref$chr, ref$start, sep=":"),
  chr = ref$chr,
  position = ref$start
)

for(s in sample_cols) {
  idx <- match(paste(methylation_dt$chr, methylation_dt$position), 
               paste(ref$chr, ref$start))
  M_col <- M_matrix[, s]
  Cov_col <- Cov_matrix[, s]
  methylation_dt[, (s) := ifelse(Cov_col > 0, M_col/Cov_col, NA_real_)]
}

coverage_dt <- data.table(
  cpg_id = paste(ref$chr, ref$start, sep=":"),
  chr = ref$chr,
  position = ref$start
)

for(s in sample_cols) {
  coverage_dt[, (s) := Cov_matrix[, s]]
}

saveRDS(methylation_dt, paste0(outdatadir, "methylation_dt.RDS"))
saveRDS(coverage_dt, paste0(outdatadir, "coverage_dt.RDS"))


cov_matrix <- getCoverage(BS)
samples_covered <- rowSums(cov_matrix >= 5, na.rm=TRUE)
keep_cpgs <- samples_covered == 4
BS_filtered <- BS[keep_cpgs]
cat("CpGs after filtering:", nrow(BS_filtered), "\n")




# ============================================================
# DMC CALLING WITH DSS
# ============================================================

group1_samples <- sample_info[condition == "normal", sample_id]
group2_samples <- sample_info[condition == "tumor", sample_id]

if(!file.exists(paste0(outdatadir, "dss_dmlTest_result.RDS"))) {
  dmlTest <- DMLtest(BS_filtered,
                     group1 = group1_samples,
                     group2 = group2_samples,
                     ncores = 6,
                     smoothing = FALSE)
  saveRDS(dmlTest, paste0(outdatadir, "dss_dmlTest_result.RDS"))
} else {
  dmlTest <- readRDS(paste0(outdatadir, "dss_dmlTest_result.RDS"))
}

# Call DMCs p <0.05
dmcs_sig <- callDML(dmlTest, p.threshold=0.05)
dmcs_sig_dt <- as.data.table(dmcs_sig)

# Call DMCs p <0.05, deltaB >0.1
dmcs_strict <- callDML(dmlTest, p.threshold=0.05, delta=0.1)
dmcs_strict_dt <- as.data.table(dmcs_strict)

dmcs_strict_dt[, direction := ifelse(diff > 0, "Hypomethylated", "Hypermethylated")]

#trecias klausimas
cat("DMCs (p<0.05):", nrow(dmcs_sig_dt), "\n")
cat("DMCs (p<0.05 + delta>0.1):", nrow(dmcs_strict_dt), "\n")

#Usage:

#     callDML(DMLresult, delta=0.1, p.threshold=1e-5)
    

#tik pasitikrinu....
dmcs_strict2 <- callDML(dmlTest, p.threshold=0.05, delta=0.2)
cat("DMCs (delta>0.2):", nrow(as.data.table(dmcs_strict2)), "\n")
dmcs_strict3 <- callDML(dmlTest, p.threshold=0.05, delta=0.3)
cat("DMCs (delta>0.3):", nrow(as.data.table(dmcs_strict3)), "\n")



# Išsaugojam
write.csv(dmcs_sig_dt, paste0(outdatadir, "dmcs_significant.csv"), row.names=FALSE)
write.csv(dmcs_strict_dt, paste0(outdatadir, "dmcs_strict.csv"), row.names=FALSE)
saveRDS(dmcs_sig_dt, paste0(outdatadir, "dmcs_sig.RDS"))
saveRDS(dmcs_strict_dt, paste0(outdatadir, "dmcs_strict.RDS"))

n_hyper <- sum(dmcs_strict_dt$diff < 0)
n_hypo <- sum(dmcs_strict_dt$diff > 0)
cat("Hypermethylated:", n_hyper, "\n")
cat("Hypomethylated:", n_hypo, "\n")


#wilcoxon
library(data.table)

normal_samples <- c("SRR11647648", "SRR11647649")
tumor_samples  <- c("SRR11647658", "SRR11647659")

methylation_dt <- readRDS(paste0(outdatadir, "methylation_dt.RDS"))
coverage_dt    <- readRDS(paste0(outdatadir, "coverage_dt.RDS"))

# Filtras
cov_filter <- coverage_dt[, SRR11647648 >= 5 & 
                             SRR11647649 >= 5 & 
                             SRR11647658 >= 5 & 
                             SRR11647659 >= 5]

meth_filtered <- methylation_dt[cov_filter]
cat("CpG po filtravimo:", nrow(meth_filtered), "\n")

# Matrica
mat <- as.matrix(meth_filtered[, c(normal_samples, tumor_samples), with=FALSE])

# Wilcoxon 
results <- row_wilcoxon_twosample(mat[, 1:2], mat[, 3:4])

# Rezultatai
wilcox_results <- data.table(
  cpg_id      = meth_filtered$cpg_id,
  chr         = meth_filtered$chr,
  position    = meth_filtered$position,
  mean_normal = rowMeans(mat[, 1:2], na.rm=TRUE),
  mean_tumor  = rowMeans(mat[, 3:4], na.rm=TRUE),
  pvalue      = results$pvalue
)

wilcox_results[, delta := mean_tumor - mean_normal]
wilcox_results[, padj  := p.adjust(pvalue, method="BH")]

# DMC calling
wilcox_dmcs <- wilcox_results[!is.na(padj) & padj < 0.05 & abs(delta) > 0.1]

cat("Wilcoxon DMCs:", nrow(wilcox_dmcs), "\n")

saveRDS(wilcox_dmcs, "/home/genetics/homework2/output/dmc_analysis/wilcox_dmcs.RDS")




# ============================================================
# COMPARE DSS VS WILCOXON
# ============================================================

# Create CpG IDs
dss_ids <- paste(dmcs_strict_dt$chr,
                 dmcs_strict_dt$pos,
                 sep = "_")

wilcox_ids <- paste(wilcox_dmcs$chr,
                    wilcox_dmcs$position,
                    sep = "_")

# If no Wilcoxon DMCs
if(length(wilcox_ids) == 0){
  wilcox_ids <- "No_DMCs"
}

# Overlap
common_dmcs <- intersect(dss_ids, wilcox_ids)

cat("Common DMCs:", length(common_dmcs), "\n")

# ============================================================
# COMPARISON TABLE
# ============================================================

comparison_table <- data.table(
  Method = c("DSS", "Wilcoxon"),

  Number_of_DMCs = c(
    nrow(dmcs_strict_dt),
    nrow(wilcox_dmcs)
  ),

  Mean_effect_size = c(
    mean(abs(dmcs_strict_dt$diff), na.rm = TRUE),
    mean(abs(wilcox_dmcs$delta), na.rm = TRUE)
  )
)

print(comparison_table)

write.csv(
  comparison_table,
  paste0(outdatadir, "method_comparison.csv"),
  row.names = FALSE
)



# ============================================================
# THRESHOLD EXPLORATION (DSS)
# ============================================================

dmcs_05_01 <- as.data.table(callDML(dmlTest, p.threshold=0.05, delta=0.1))
dmcs_01_01 <- as.data.table(callDML(dmlTest, p.threshold=0.01, delta=0.1))
dmcs_05_02 <- as.data.table(callDML(dmlTest, p.threshold=0.05, delta=0.2))
dmcs_01_02 <- as.data.table(callDML(dmlTest, p.threshold=0.01, delta=0.2))
all_dml <- as.data.table(dmlTest)

dmcs_05_01 <- all_dml[pval < 0.05 & abs(diff) > 0.1]
dmcs_01_01 <- all_dml[pval < 0.01 & abs(diff) > 0.1]
dmcs_05_02 <- all_dml[pval < 0.05 & abs(diff) > 0.2]
dmcs_01_02 <- all_dml[pval < 0.01 & abs(diff) > 0.2]

threshold_table <- data.table(
  Threshold = c(
    "p<0.05 & |Δβ|>0.1",
    "p<0.01 & |Δβ|>0.1",
    "p<0.05 & |Δβ|>0.2",
    "p<0.01 & |Δβ|>0.2"
  ),

  Number_of_DMCs = c(
    nrow(dmcs_05_01),
    nrow(dmcs_01_01),
    nrow(dmcs_05_02),
    nrow(dmcs_01_02)
  )
)

threshold_table


ggplot(threshold_table, aes(x=Threshold, y=Number_of_DMCs, group=1)) +
  geom_line(color="#332932", linewidth=1) +
  geom_point(color="#e89eb8", size=4) +
  geom_text(aes(label=Number_of_DMCs), vjust=-1, size=4) +
  labs(title="Number of DMCs vs Threshold Stringency",
       x="Threshold", y="Number of DMCs") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(paste0(outdatadir, "dmc_threshold_plot.png"), width=10, height=6)


# ============================================================
# ANNOTATION
# ============================================================

if(!file.exists(paste0(outdatadir, "annotation_hg19.RDS"))) {
  annots <- c('hg19_cpgs',
              'hg19_genes_promoters',
              'hg19_genes_5UTRs',
              'hg19_genes_exons',
              'hg19_genes_introns',
              'hg19_genes_3UTRs',
              'hg19_genes_intergenic')
  annotations <- build_annotations(genome='hg19', annotations=annots)
  saveRDS(annotations, paste0(outdatadir, "annotation_hg19.RDS"))
} else {
  annotations <- readRDS(paste0(outdatadir, "annotation_hg19.RDS"))
}

# Convert DMCs to GRanges - chr prefix needed for hg19
dmcs_gr <- GRanges(
  seqnames = paste0("chr", dmcs_strict_dt$chr),
  ranges = IRanges(start=dmcs_strict_dt$pos, width=1),
  diff = dmcs_strict_dt$diff,
  pval = dmcs_strict_dt$pval,
  fdr = dmcs_strict_dt$fdr,
  direction = dmcs_strict_dt$direction
)

# Annotate
dmcs_annotated <- annotate_regions(
  regions = dmcs_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE
)

dmcs_annot_dt <- as.data.table(dmcs_annotated)
dmcs_annot_dt[, annot_type := gsub("hg19_", "", annot.type)]
dmcs_annot_dt[, annot_type := gsub("genes_", "", annot_type)]
dmcs_annot_dt[, annot_type := gsub("cpg_", "", annot_type)]

# Promoter DMCs
promoter_dmcs <- dmcs_annot_dt[grepl("promoter", annot_type, ignore.case=TRUE)]

cat("Total DMCs:", nrow(dmcs_strict_dt), "\n")
cat("Promoter DMCs:", nrow(promoter_dmcs), "\n")
cat("% in promoters:", round(nrow(promoter_dmcs)/nrow(dmcs_strict_dt)*100, 1), "%\n")

write.csv(promoter_dmcs, paste0(outdatadir, "promoter_dmcs.csv"), row.names=FALSE)
saveRDS(dmcs_annot_dt, paste0(outdatadir, "dmcs_annotated.RDS"))

# ============================================================
# VOLCANO PLOT
# ============================================================

volcano_df <- as.data.table(dmlTest)
volcano_df[, significant := ifelse(pval < 0.05 & abs(diff) > 0.1, "Significant", "Not significant")]

ggplot(volcano_df, aes(x=diff, y=-log10(pval), color=significant)) +
  geom_point(alpha=0.5, size=0.8) +
  scale_color_manual(values=c("Not significant"="grey", "Significant"="#332932")) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="red") +
  geom_vline(xintercept=c(-0.1, 0.1), linetype="dashed", color="blue") +
  labs(title="Volcano Plot of DMCs",
       x="Methylation Difference (Δβ)",
       y="-log10(p-value)") +
  theme_bw()

ggsave(paste0(outdatadir, "volcano_plot.png"), width=10, height=7)

# ============================================================
# MA PLOT
# ============================================================
meth_means <- rowMeans(
  methylation_dt[, c("SRR11647648","SRR11647649","SRR11647658","SRR11647659"), with=FALSE],
  na.rm=TRUE)

volcano_df[, mean_methylation := meth_means[match(
  paste(volcano_df$chr, volcano_df$pos),
  paste(methylation_dt$chr, methylation_dt$position)
)]]

ggplot(volcano_df, aes(x=mean_methylation, y=diff, color=significant)) +
  geom_point(alpha=0.3, size=0.5) +
  scale_color_manual(values=c("Not significant"="grey80", "Significant"="#332932")) +
  geom_hline(yintercept=0, linetype="dashed", color="red") +
  labs(title="MA Plot of DMCs",
       x="Average Methylation Level",
       y="Methylation Difference (Δβ)") +
  theme_bw()

ggsave(paste0(outdatadir, "ma_plot.png"), width=10, height=7)
#




#sestas sklausimas
cat("Total DMCs:", nrow(dmcs_strict_dt), "\n")
cat("Promoter DMCs:", nrow(promoter_dmcs), "\n")
cat("% in promoters:", round(nrow(promoter_dmcs)/nrow(dmcs_strict_dt)*100, 1), "%\n")











# ============================================================
# DMR CALLING
# ============================================================

# Parameter set 1
dmr1 <- callDMR(dmlTest, 
                p.threshold=0.05,
                minlen=50,
                minCG=3,
                dis.merge=100)
dmr1_dt <- as.data.table(dmr1)
dmr1_dt[, cpg_density := nCG / length]

# Parameter set 2
dmr2 <- callDMR(dmlTest,
                p.threshold=0.01,
                minlen=100,
                minCG=5,
                dis.merge=100)
dmr2_dt <- as.data.table(dmr2)
dmr2_dt[, cpg_density := nCG / length]

# Parameter set 3
dmr3 <- callDMR(dmlTest,
                p.threshold=0.01,
                minlen=200,
                minCG=7,
                dis.merge=100)
dmr3_dt <- as.data.table(dmr3)
dmr3_dt[, cpg_density := nCG / length]

# Comparison table
dmr_comparison <- data.table(
  Parameters = c("minCG=3, minlen=50, p=0.05",
                 "minCG=5, minlen=100, p=0.01",
                 "minCG=7, minlen=200, p=0.01"),
  N_DMRs = c(nrow(dmr1_dt), nrow(dmr2_dt), nrow(dmr3_dt)),
  Mean_length = c(mean(dmr1_dt$length), mean(dmr2_dt$length), mean(dmr3_dt$length)),
  Mean_CpG_density = c(mean(dmr1_dt$cpg_density), mean(dmr2_dt$cpg_density), mean(dmr3_dt$cpg_density))
)

print(dmr_comparison)
write.csv(dmr_comparison, paste0(outdatadir, "dmr_comparison.csv"), row.names=FALSE)

# Save all DMR sets
saveRDS(dmr1_dt, paste0(outdatadir, "dmr1.RDS"))
saveRDS(dmr2_dt, paste0(outdatadir, "dmr2.RDS"))
saveRDS(dmr3_dt, paste0(outdatadir, "dmr3.RDS"))


# ============================================================
# DMR CHARACTERIZATION
# ============================================================

dmr <- dmr1_dt
dmr[, direction := ifelse(diff.Methy > 0, "Hypomethylated", "Hypermethylated")]

# 1. Length distribution
p1 <- ggplot(dmr, aes(x=length, fill=direction)) +
  geom_histogram(bins=50, alpha=0.7) +
  scale_fill_manual(values=c("Hypomethylated"="#e89eb8", "Hypermethylated"="#332932")) +
  labs(title="DMR Length Distribution", x="Length (bp)", y="Count") +
  theme_bw()

# 2. CpGs per DMR
p2 <- ggplot(dmr, aes(x=nCG, fill=direction)) +
  geom_histogram(bins=30, alpha=0.7) +
  scale_fill_manual(values=c("Hypomethylated"="#e89eb8", "Hypermethylated"="#332932")) +
  labs(title="CpGs per DMR", x="Number of CpGs", y="Count") +
  theme_bw()

# 3. Hyper vs Hypo
p3 <- ggplot(dmr, aes(x=direction, fill=direction)) +
  geom_bar() +
  geom_text(stat="count", aes(label=after_stat(count)), vjust=-0.5) +
  scale_fill_manual(values=c("Hypomethylated"="#e89eb8", "Hypermethylated"="#332932")) +
  labs(title="Hyper vs Hypo DMRs", x="Direction", y="Count") +
  theme_bw()

# 4. Density plot across genome
p4 <- ggplot(dmr, aes(x=start, fill=direction)) +
  geom_density(alpha=0.5) +
  scale_fill_manual(values=c("Hypomethylated"="#e89eb8", "Hypermethylated"="#332932")) +
  labs(title="DMR Density across Chromosome 20", x="Position", y="Density") +
  theme_bw()

# 5. Correlation between DMR length and effect size
p5 <- ggplot(dmr, aes(x=length, y=abs(diff.Methy), color=direction)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm") +
  scale_color_manual(values=c("Hypomethylated"="#e89eb8", "Hypermethylated"="#332932")) +
  labs(title="DMR Length vs Effect Size",
       x="Length (bp)", y="|Mean Methylation Difference|") +
  theme_bw()

# Combine all plots
library(patchwork)
combined <- (p1 + p2) / (p3 + p4) / p5
ggsave(paste0(outdatadir, "dmr_characterization.png"), combined, width=14, height=15)



# Save DMR list
dmr_final <- dmr[, .(
  chr,
  start,
  end,
  length,
  nCG,
  meanMethy1,
  meanMethy2,
  diff.Methy,
  areaStat,
  direction
)]

write.csv(dmr_final, paste0(outdatadir, "dmr_final_list.csv"), row.names=FALSE)
saveRDS(dmr_final, paste0(outdatadir, "dmr_final.RDS"))

cat("DMR list saved!\n")
cat("Total DMRs:", nrow(dmr_final), "\n")



# ============================================================
# COMPARE DMCs AND DMRs
# ============================================================

library(GenomicRanges)

# Convert DMCs and DMRs to GRanges
dmc_gr <- GRanges(
  seqnames = paste0("chr", dmcs_strict_dt$chr),
  ranges = IRanges(start=dmcs_strict_dt$pos, width=1)
)

dmr_gr <- GRanges(
  seqnames = paste0("chr", dmr$chr),
  ranges = IRanges(start=dmr$start, end=dmr$end)
)

# DMCs falling within DMRs
overlaps <- findOverlaps(dmc_gr, dmr_gr)
dmcs_in_dmr <- length(unique(queryHits(overlaps)))
dmcs_isolated <- nrow(dmcs_strict_dt) - dmcs_in_dmr

cat("DMCs in DMRs:", dmcs_in_dmr, "\n")
cat("Isolated DMCs:", dmcs_isolated, "\n")
cat("% DMCs in DMRs:", round(dmcs_in_dmr/nrow(dmcs_strict_dt)*100, 1), "%\n")

# Isolated DMRs (no DMCs inside - shouldn't happen but check)
dmrs_with_dmcs <- length(unique(subjectHits(overlaps)))
isolated_dmrs <- nrow(dmr) - dmrs_with_dmcs
cat("Isolated DMRs:", isolated_dmrs, "\n")

# Plot DMC fractions
frac_df <- data.frame(
  category = c("DMCs in DMRs", "Isolated DMCs"),
  count = c(dmcs_in_dmr, dmcs_isolated)
)

# Plot enrichment
enrichment_df <- data.frame(
  category = c("DMCs in DMRs", "Isolated DMCs", 
                "DMRs in Promoters", "DMRs in CGI", "Isolated DMRs"),
  count = c(dmcs_in_dmr, dmcs_isolated,
            length(promoter_dmrs), length(cgi_dmrs), isolated_dmrs),
  type = c("DMC", "DMC", "DMR", "DMR", "DMR")
)

ggplot(enrichment_df, aes(x=category, y=count, fill=type)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=count), vjust=-0.5, size=4) +
  scale_fill_manual(values=c("DMC"="#e89eb8", "DMR"="#332932")) +
  labs(title="DMC and DMR Overlap Summary",
       x="Category", y="Count") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(paste0(outdatadir, "dmc_dmr_summary.png"), width=10, height=6)
# Enrichment in promoters and CGI
promoter_dmrs <- subsetByOverlaps(dmr_gr, annotations[grepl("promoter", annotations$type)])
cgi_dmrs <- subsetByOverlaps(dmr_gr, annotations[grepl("island", annotations$type)])

cat("DMRs in promoters:", length(promoter_dmrs), "\n")
cat("DMRs in CGI:", length(cgi_dmrs), "\n")
cat("% DMRs in promoters:", round(length(promoter_dmrs)/nrow(dmr)*100, 1), "%\n")
cat("% DMRs in CGI:", round(length(cgi_dmrs)/nrow(dmr)*100, 1), "%\n")


# ============================================================
# HEATMAP OF TOP DMCs AND DMRs
# ============================================================

library(pheatmap)

# Top 100 DMCs by p-value
top_dmcs <- dmcs_strict_dt[order(pval)][1:100]
top_dmcs[, cpg_id := paste(chr, pos, sep=":")]

# Get methylation values for top DMCs
meth_top <- methylation_dt[cpg_id %in% top_dmcs$cpg_id]
meth_matrix <- as.matrix(meth_top[, c("SRR11647648","SRR11647649","SRR11647658","SRR11647659"), with=FALSE])
rownames(meth_matrix) <- meth_top$cpg_id

# Remove NA rows
meth_matrix <- meth_matrix[complete.cases(meth_matrix),]

# Annotation for columns
annotation_col <- data.frame(
  Condition = c("Normal", "Normal", "Tumor", "Tumor"),
  row.names = c("SRR11647648","SRR11647649","SRR11647658","SRR11647659")
)

ann_colors <- list(
  Condition = c(Normal="#e89eb8", Tumor="#332932")
)

# Plot
png(paste0(outdatadir, "heatmap_top_dmcs.png"), width=1000, height=1200)
pheatmap(meth_matrix,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("blue","white","red"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         main = "Top 100 DMCs by p-value")
dev.off()
