
#BiocManager::install("ChIPseeker")
library("ChIPseeker")

#BiocManager::install("GenomicFeatures")
library("GenomicFeatures")

#install.packages("ggplot2")
library("ggplot2")

#BiocManager::install("txdbmaker")
library("txdbmaker")

setwd("/home/gyk/project/lw_atac_nf/")

file_annot  <- "/home/gyk/reference/gff/Bomo_gene_models_chr_update.gff3"
file_format <- "gff"
file_peaks  <- "results/bwa/merged_replicate/macs2/narrow_peak/DZ.mRp.clN_peaks.narrowPeak"
sample_name <- "DZ"
output_dir  <- "analysis/results/peaks_annotation/"
flag_chrs   <- "no"


# 使用makeTxDbFromGFF函数解析GFF或GTF文件制作TxDb
txdb       <- makeTxDbFromGFF(file = file_annot, format = file_format)

# 使用readPeakFile函数读取peaks文件，生成一个对象
data_peaks <- readPeakFile(file_peaks)

# getPromoters获取启动子上下游3kb的位置
promoter   <- getPromoters(TxDb = txdb, upstream=3000, downstream=3000)

# getTagMatrix peaks结果转换成矩阵 窗口为启动子
tagMatrix  <- getTagMatrix(data_peaks, windows = promoter)
peaksAnno  <- annotatePeak(data_peaks, TxDb=txdb, tssRegion=c(-3000, 3000), 
                           verbose=FALSE, addFlankGeneInfo=TRUE, flankDistance=5000)

# 转为表格
annot_tab <- as.data.frame(peaksAnno)
colnames(annot_tab)[1] <- "#seqnames"

write.table(x = annot_tab, file = paste0(output_dir, "/", sample_name, ".peaks.gene_annotation.xls"),
            quote = FALSE, sep = "\t", row.names = FALSE)
#cat(paste0(paste(unique(annot_tab$geneId),collapse="\n"),"\n"), file=paste0(output_dir, "/", sample_name, ".peaks.gene_annotation.list"))

# 绘图 
# 可视化peak文件 横坐标:染色体位点; 纵坐标:每条染色体上的peaks位点及得分
pdf(file = paste0(output_dir, "/", sample_name, ".peaks.chromosome_distribution.pdf"))
covplot(data_peaks, weightCol = 'V5', xlab = "Chromosome Size (bp)", 
        title = paste0(sample_name, " Peaks over Chromosomes"))
dev.off()

# tagHeatmap 函数查看TSS(转录起始位点)上下游3kb区域peak的分布情况（热图）
pdf(file = paste0(output_dir, "/", sample_name, ".peaks.around_tss_heatmap.pdf"))
tagHeatmap(tagMatrix, xlab="Genomic Region (5'->3')", title = paste0(sample_name," peak distribution arround TSS"))
dev.off()

# plotAvgProf 函数绘制一个峰图描述所有分布的平均情况
pdf(file = paste0(output_dir, "/", sample_name, ".peaks.around_tss_profile.pdf"), 6, 4)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", 
            ylab = "Peak Count Frequency", title = paste0(sample_name," peak distribution arround TSS"))
dev.off()

# plotDistToTSS 函数计算出距离最近的基因的TSS上游和下游结合位点的百分比,并可视化其分布
pdf(file = paste0(output_dir, "/", sample_name, ".peaks.around_tss_bar.pdf"), 7, 3)
plotDistToTSS(peaksAnno,title="Distribution of transcription factor-binding loci \n relative to TSS")
dev.off()

# plotAnnoBar 函数可视化peaks组成比例:柱状图
pdf(file = paste0(output_dir, "/", sample_name, ".peaks.feature_distribution_bar.pdf"), 7, 3)
plotAnnoBar(x = peaksAnno, ylab = "Percentage(%)", title = paste0(sample_name, " Feature Distribution"))
dev.off()

# plotAnnoPie 函数可视化peaks组成比例:饼图
pdf(file = paste0(output_dir, "/", sample_name, ".peaks.feature_distribution_pie.pdf"), 6, 4)
plotAnnoPie(x = peaksAnno, pie3D = FALSE)
dev.off()

#__END__
