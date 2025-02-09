# BiocManager::install("DESeq2") if necessary

# Set the working directory to the directory containing the data
setwd("/home/gyk/project/lw_diffanalysis_nf/")

DEG <- read.table("results/tables/differential/condition_control_treated.deseq2.results.tsv", header = TRUE)
DEG <- na.omit(DEG)

# Volcano plot
logFC_cutoff <- with(DEG, mean(abs(log2FoldChange)) + 2 * sd(abs(log2FoldChange)))
# logFC_cutoff  <- 1.5
DEG$change <- as.factor(ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                               ifelse(DEG$log2FoldChange > logFC_cutoff, "UP", "DOWN"), "NOT")
)
write.csv(DEG, "analysis/results/DEG_deseq2.results.csv")
this_tile <- paste0("Cutoff for logFc is ", round(logFC_cutoff, 3),
  "\nThe number of up-regulated genes is ", nrow(DEG[DEG$change == "UP", ]),
  "\nThe number of down-regulated genes is ", nrow(DEG[DEG$change == "DOWN", ])
)

library(ggplot2)
g <- ggplot(data = DEG,
            aes(x = log2FoldChange, y = -log10(padj),
                color = change)) +
  geom_point(alpha = 0.4, size = 1.75) +
  theme_set(theme_set(theme_bw(base_size = 20))) +
  xlab("log2 fold change") + ylab("-log10 padj") +
  ggtitle(this_tile) + theme(plot.title = element_text(size = 15, hjust = 0.5)) +
  scale_colour_manual(values = c("blue", "black", "red")) +
  # corresponding to the levels(res$change)
  geom_hline(yintercept = -log10(0.05), lty = 4) + # Define p-value and line type
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 4) # Define fold change and line type
print(g)
ggsave(g, filename = "analysis/results/volcano.pdf", device = "pdf", dpi = 720, width = 6, height = 7)


# pheatmap
library(dplyr)
geneid <- filter(DEG, change != "NOT") %>% arrange(padj) %>% select(gene_id)
write.table(geneid, "analysis/results/DEG_id.list", row.names = FALSE, col.names = FALSE, quote = FALSE)
tpm <- read.table("/home/gyk/project/lw_rna_nf/results/star_salmon/salmon.merged.gene_tpm.tsv", header = TRUE)
nrDEG <- tpm %>%
  filter(gene_id %in% geneid$gene_id) %>%
  mutate(gene_id = factor(gene_id, levels = geneid$gene_id)) %>%
  arrange(gene_id) %>%
  select(DZ_rep1, DZ_rep2, DZ_rep3, Ava_rep1, Ava_rep2, Ava_rep3)
rownames(nrDEG) <- geneid$gene_id
# install.packages("pheatmap") if necessary
library(pheatmap)
choose_gene <- head(rownames(nrDEG), 50)##50 maybe better
choose_matrix <- nrDEG[choose_gene, ]

# load annotation file
annotation <- read.table("analysis/Bmo_annotations.tsv", sep = "\t", header = TRUE, row.names = 1)
choose_matrix_annotation <- annotation[choose_gene, ]
rownames(choose_matrix_annotation) <- choose_gene
# replace NA to -
choose_matrix_annotation[is.na(choose_matrix_annotation)] <- "-"

# Log-transform the choose_matrix matrix
log_choose_matrix <- log2(choose_matrix + 1)  # Add 1 to avoid taking the logarithm of 0

# plot heatmap for top 50 genes on padj
pheatmap(
  log_choose_matrix,
  # annotation_col = choose_matrix_annotation["Gene_Name"],
  annotation_row = choose_matrix_annotation["Gene_Name"],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 8,
  cellwidth = 20,
  cellheight = 10,
  main = "Gene Expression Heatmap",
  legend = TRUE,
  filename = "analysis/results/DEG_top50_heatmap.pdf"
)

# plot heatmap for all genes of nrDEG
log_nrDEG <- log2(nrDEG + 1)
pheatmap(
  log_nrDEG,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize = 8,
  cellwidth = 30,
  cellheight = 1,
  main = "Gene Expression Heatmap",
  legend = TRUE,
  filename = "analysis/results/DEG_heatmap.pdf"
)
