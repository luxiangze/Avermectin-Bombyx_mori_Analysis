#!/usr/bin/env Rscript

# usage: Rscript enrichx.R -i gene_list.txt -o ./enrichment_results -d org.Bmori.eg.db

# Parse command line arguments
library(optparse, quietly = TRUE)

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input gene list file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="./enrichment_results", 
              help="Output directory [default= %default]", metavar="character"),
  make_option(c("-d", "--database"), type="character", default="org.Bmori.eg.db", 
              help="Organism database [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file must be supplied", call.=FALSE)
}

# Get input file prefix for output files
file_prefix <- tools::file_path_sans_ext(basename(opt$input))

# R包安装
# if(!require(BiocManager)){
# 	install.packages("BiocManager");
# }
# if(!require(org.Hs.eg.db)){
# 	BiocManager::install("org.Hs.eg.db")
# }
# library(org.Hs.eg.db)

if(!require(org.Bmori.eg.db)){
  install.packages("analysis/org.Bmori.eg.db_1.0.tar.gz")
}
library(org.Bmori.eg.db, quietly = TRUE)

# if(!require(clusterProfiler)){
# 	BiocManager::install("clusterProfiler")
# }
library(clusterProfiler, quietly = TRUE)

# if(!require(ggplot2)){
# 	install.packages("ggplot2")
# }
library(ggplot2, quietly = TRUE)

# if(!require(dplyr)){
#   install.packages("dplyr")
# }
library(dplyr, quietly = TRUE)

# Create output directory if it doesn't exist
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Read input gene list
gene_list <- readLines(opt$input)

# GO富集分析
go_enrich_result <- enrichGO(gene = gene_list, 
    keyType = "GID", 
    ont = "ALL",               # Sub Ontology
    OrgDb = get(opt$database), # 物种数据库
    pAdjustMethod = "fdr",     # p值的校正方法
    pvalueCutoff = 1,       # p值的阈值
    qvalueCutoff = 1)       # q值的阈值

# 将GO富集结果(数据框)转换为表格
go_table <- as.data.frame(go_enrich_result)

# 计算富集因子
split_col4 <- strsplit(go_table$GeneRatio, "/")
split_col5 <- strsplit(go_table$BgRatio, "/")
number_col4 <- sapply(split_col4, function(x) as.numeric(x))
number_col5 <- sapply(split_col5, function(x) as.numeric(x))
Enrich_factor <- number_col4[1, ] * number_col5[2, ] / number_col4[2, ] / number_col5[1, ]
GO_table <- cbind(go_table[, 1:5], Enrich_factor=Enrich_factor, go_table[, 6:ncol(go_table)])

# 输出GO分析结果
write.table(GO_table, 
            file=file.path(opt$outdir, paste0(file_prefix, "_go_enrichment.xls")), 
            sep="\t", quote=FALSE, row.names=FALSE)

# GO柱状图
pdf(file=file.path(opt$outdir, paste0(file_prefix, "_go_barplot.pdf")), width=13, height=10) 
barplot(go_enrich_result, showCategory=5, split="ONTOLOGY", col="p.adjust", 
        label_format=35, title="GO Enrichment Barplot", font.size=15) +
       facet_grid(ONTOLOGY~., scale='free')
dev.off()

# GO气泡图
pdf(file=file.path(opt$outdir, paste0(file_prefix, "_go_bubble.pdf")), width=13, height=10) 
dotplot(go_enrich_result, showCategory=5, split="ONTOLOGY", col="p.adjust", 
        label_format=35, title="GO Enrichment Bubble", font.size=15) +
        facet_grid(ONTOLOGY~., scale='free')
dev.off()

# GO富集因子图
BP <- subset(GO_table, ONTOLOGY=="BP")[1:5, ]
CC <- subset(GO_table, ONTOLOGY=="CC")[1:5, ]
MF <- subset(GO_table, ONTOLOGY=="MF")[1:5, ]
topgo <- rbind(BP, CC, MF)

pdf(file=file.path(opt$outdir, paste0(file_prefix, "_go_enrichfactor.pdf")), width=13, height=10)
ggplot(topgo, aes(Enrich_factor, Description)) + 
        geom_point(aes(size=Count, color=p.adjust)) +
        scale_color_gradient(low="red", high="blue") + 
        labs(color="p.adjust", size="Count", x="Enrich_factor", y="GO term", 
             title="Bubble of GO Enrich Factor") +
        theme_bw() +
        theme(
              plot.title  = element_text(size=15),
              axis.text.x = element_text(size=15),
              axis.text.y = element_text(size=15),
              axis.title.x = element_text(size=15),
              axis.title.y = element_text(size=15)
        ) +
        facet_grid(ONTOLOGY~., scale='free')
dev.off()

# KEGG富集分析
input_data <- read.table(opt$input, header = FALSE, stringsAsFactors = FALSE)
all_enretz_ids <- read.table("analysis/entrez.txt", header = TRUE, sep = "\t")
entrez_ids <- all_enretz_ids %>% filter(all_enretz_ids$GENE %in% input_data$V1) %>% select(ENTREZ.ID)
entrez_ids <- as.character(entrez_ids$ENTREZ.ID)

kegg_enrich_result <- enrichKEGG(
    gene = entrez_ids,
    keyType = "kegg",
    organism = "bmor",  # 物种
    pAdjustMethod = "fdr",  # p值校正方法
    pvalueCutoff = 1,
    qvalueCutoff = 1)
    
# 将KEGG富集结果转换为表格
kegg_table <- as.data.frame(kegg_enrich_result)
split_col4 <- strsplit(kegg_table$GeneRatio, "/")
split_col5 <- strsplit(kegg_table$BgRatio, "/")
number_col4 <- sapply(split_col4, function(x) as.numeric(x))
number_col5 <- sapply(split_col5, function(x) as.numeric(x))
Enrich_factor <- number_col4[1, ] * number_col5[2, ] / number_col4[2, ] / number_col5[1, ]
KEGG_table <- cbind(kegg_table[, 1:5], Enrich_factor=Enrich_factor, kegg_table[, 6:ncol(kegg_table)])

# 输出KEGG分析结果
write.table(KEGG_table, 
            file=file.path(opt$outdir, paste0(file_prefix, "_kegg_enrichment.xls")), 
            sep="\t", quote=FALSE, row.names=FALSE)

# KEGG柱状图
pdf(file=file.path(opt$outdir, paste0(file_prefix, "_kegg_barplot.pdf")), width=13, height=10)
barplot(kegg_enrich_result, showCategory=15, title="KEGG Enrichment Barplot", 
        col="p.adjust", label_format=35, font.size=15)
dev.off()

# KEGG气泡图
pdf(file=file.path(opt$outdir, paste0(file_prefix, "_kegg_bubble.pdf")), width=13, height=10)
dotplot(kegg_enrich_result, showCategory=15, title="KEGG Enrichment Bubble", 
        col="p.adjust", label_format=35, font.size=15)
dev.off()

# KEGG富集因子图
topkegg <- head(KEGG_table, 15)
# 删除description中" - ""后面全部的内容
topkegg <- topkegg %>% mutate(Description = gsub(" - .*", "", Description))
pdf(file=file.path(opt$outdir, paste0(file_prefix, "_kegg_enrichfactor.pdf")), width=13, height=10)
ggplot(topkegg, aes(Enrich_factor, Description)) + 
        geom_point(aes(size=Count, color=p.adjust)) +
        scale_color_gradient(low="red", high="blue") + 
        labs(color="p.adjust", size="Count", x="Enrich_factor", y="KEGG Pathway", 
             title="Bubble of KEGG Enrich Factor") +
        theme_bw() +
        theme(
              plot.title = element_text(size=15),
              axis.text.x = element_text(size=15),
              axis.text.y = element_text(size=15),
              axis.title.x = element_text(size=15),
              axis.title.y = element_text(size=15)
        )
dev.off()

# __END__
