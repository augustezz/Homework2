# ============================================================
# enrichment.R
# ============================================================

library(pacman)
p_load(data.table, clusterProfiler, org.Hs.eg.db, ggplot2)

# Direktorijos
indatadir <- "/home/genetics/homework2/output/dmc_analysis/"
outdatadir <- "/home/genetics/homework2/output/enrich/"
dir.create(outdatadir, recursive=TRUE)

# Įkelk duomenis
dmcs_annot_dt <- readRDS(paste0(indatadir, "dmcs_annotated.RDS"))
dmcs_annot_dt <- as.data.table(dmcs_annot_dt)

# Patikrink kolonas
cat("Kolonos:", colnames(dmcs_annot_dt), "\n")

# Tvarkome kolonų pavadinimus
dmcs_annot_dt[, annot_type := gsub("hg19_", "", annot.type)]
dmcs_annot_dt[, annot_type := gsub("genes_", "", annot_type)]

# Patikrink ar direction yra
if(!"direction" %in% colnames(dmcs_annot_dt)) {
  dmcs_annot_dt[, direction := ifelse(diff > 0, "Hypomethylated", "Hypermethylated")]
}

# Genų sąrašai
all_dmc_genes     <- unique(dmcs_annot_dt[!is.na(annot.gene_id), annot.gene_id])
promoter_dmc_genes <- unique(dmcs_annot_dt[grepl("promoter", annot_type) & 
                                            !is.na(annot.gene_id), annot.gene_id])
hyper_genes <- unique(dmcs_annot_dt[direction == "Hypermethylated" & 
                                     !is.na(annot.gene_id), annot.gene_id])
hypo_genes  <- unique(dmcs_annot_dt[direction == "Hypomethylated" & 
                                     !is.na(annot.gene_id), annot.gene_id])

cat("All DMC genes:", length(all_dmc_genes), "\n")
cat("Promoter DMC genes:", length(promoter_dmc_genes), "\n")
cat("Hyper genes:", length(hyper_genes), "\n")
cat("Hypo genes:", length(hypo_genes), "\n")

# GO enrichment funkcija
run_go <- function(genes, label) {
  cat("Running GO for:", label, "\n")
  if(length(genes) < 5) { cat(label, "- per mažai genų\n"); return(NULL) }
  tryCatch(
    enrichGO(gene = genes,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2,
             readable = TRUE),
    error = function(e) { cat("Klaida:", e$message, "\n"); return(NULL) }
  )
}

# KEGG enrichment funkcija
run_kegg <- function(genes, label) {
  cat("Running KEGG for:", label, "\n")
  if(length(genes) < 5) { cat(label, "- per mažai genų\n"); return(NULL) }
  tryCatch(
    enrichKEGG(gene = genes,
               organism = "hsa",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2),
    error = function(e) { cat("Klaida:", e$message, "\n"); return(NULL) }
  )
}

# Paleisti visas analizes
go_all      <- run_go(all_dmc_genes, "All DMC")
go_promoter <- run_go(promoter_dmc_genes, "Promoter DMC")
go_hyper    <- run_go(hyper_genes, "Hyper")
go_hypo     <- run_go(hypo_genes, "Hypo")

kegg_all      <- run_kegg(all_dmc_genes, "All DMC")
kegg_promoter <- run_kegg(promoter_dmc_genes, "Promoter DMC")
kegg_hyper    <- run_kegg(hyper_genes, "Hyper")
kegg_hypo     <- run_kegg(hypo_genes, "Hypo")

# Išsaugoti
saveRDS(go_all,        paste0(outdatadir, "go_all.RDS"))
saveRDS(go_promoter,   paste0(outdatadir, "go_promoter.RDS"))
saveRDS(go_hyper,      paste0(outdatadir, "go_hyper.RDS"))
saveRDS(go_hypo,       paste0(outdatadir, "go_hypo.RDS"))
saveRDS(kegg_all,      paste0(outdatadir, "kegg_all.RDS"))
saveRDS(kegg_promoter, paste0(outdatadir, "kegg_promoter.RDS"))
saveRDS(kegg_hyper,    paste0(outdatadir, "kegg_hyper.RDS"))
saveRDS(kegg_hypo,     paste0(outdatadir, "kegg_hypo.RDS"))

# Grafikai
plot_go <- function(go_result, title, filename) {
  if(is.null(go_result)) return(NULL)
  if(nrow(go_result) == 0) return(NULL)
  p <- dotplot(go_result, showCategory=15, title=title)
  ggsave(paste0(outdatadir, filename), p, width=10, height=8)
}

plot_go(go_all,      "GO - All DMC genes",            "go_all.png")
plot_go(go_promoter, "GO - Promoter DMC genes",        "go_promoter.png")
plot_go(go_hyper,    "GO - Hypermethylated genes",     "go_hyper.png")
plot_go(go_hypo,     "GO - Hypomethylated genes",      "go_hypo.png")

cat("Viskas baigta!\n")


# ============================================================
# DMR ENRICHMENT
# ============================================================

library(GenomicRanges)
library(annotatr)

# Įkelk DMR ir anotacijas
dmr_final <- readRDS(paste0(indatadir, "dmr_final.RDS"))
annotations <- readRDS(paste0(indatadir, "annotation_hg19.RDS"))

# Konvertuok DMR į GRanges
dmr_gr <- GRanges(
  seqnames = paste0("chr", dmr_final$chr),
  ranges = IRanges(start=dmr_final$start, end=dmr_final$end),
  diff = dmr_final$diff.Methy,
  direction = dmr_final$direction
)

# Anotuok DMR
dmr_annotated <- annotate_regions(
  regions = dmr_gr,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE
)

dmr_annot_dt <- as.data.table(dmr_annotated)
dmr_annot_dt[, annot_type := gsub("hg19_", "", annot.type)]
dmr_annot_dt[, annot_type := gsub("genes_", "", annot_type)]

# Genų sąrašai
all_dmr_genes <- unique(dmr_annot_dt[!is.na(annot.gene_id), annot.gene_id])
promoter_dmr_genes <- unique(dmr_annot_dt[grepl("promoter", annot_type) & 
                                           !is.na(annot.gene_id), annot.gene_id])

cat("All DMR genes:", length(all_dmr_genes), "\n")
cat("Promoter DMR genes:", length(promoter_dmr_genes), "\n")

# GO ir KEGG
go_dmr_all      <- run_go(all_dmr_genes, "All DMR")
go_dmr_promoter <- run_go(promoter_dmr_genes, "Promoter DMR")

kegg_dmr_all      <- run_kegg(all_dmr_genes, "All DMR")
kegg_dmr_promoter <- run_kegg(promoter_dmr_genes, "Promoter DMR")

# Išsaugoti
saveRDS(go_dmr_all,      paste0(outdatadir, "go_dmr_all.RDS"))
saveRDS(go_dmr_promoter, paste0(outdatadir, "go_dmr_promoter.RDS"))
saveRDS(kegg_dmr_all,    paste0(outdatadir, "kegg_dmr_all.RDS"))
saveRDS(kegg_dmr_promoter, paste0(outdatadir, "kegg_dmr_promoter.RDS"))

# Grafikai
plot_go(go_dmr_all,      "GO - All DMR genes",        "go_dmr_all.png")
plot_go(go_dmr_promoter, "GO - Promoter DMR genes",   "go_dmr_promoter.png")

cat("DMR enrichment baigta!\n")


# ============================================================
# TOP TERMINAI PALYGINIMUI
# ============================================================

print_top <- function(go_result, label) {
  cat("\n===", label, "===\n")
  if(is.null(go_result) || nrow(go_result) == 0) {
    cat("Nėra reikšmingų terminų\n")
    return()
  }
  print(head(go_result@result[, c("Description", "pvalue", "Count")], 15))
}

print_top(go_all,        "All DMC genes")
print_top(go_promoter,   "Promoter DMC genes")
print_top(go_hyper,      "Hypermethylated genes")
print_top(go_hypo,       "Hypomethylated genes")
print_top(go_dmr_all,    "All DMR genes")
print_top(go_dmr_promoter, "Promoter DMR genes")



cat("All DMR genes:", length(all_dmr_genes), "\n")
cat("Promoter DMR genes:", length(promoter_dmr_genes), "\n")
cat('All DMC genes:', length(all_dmc_genes), '\n')
cat("Promoter DMR genes:", length(promoter_dmc_genes), "\n")
