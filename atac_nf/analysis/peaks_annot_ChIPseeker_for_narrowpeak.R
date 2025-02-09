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


# Use the makeTxDbFromGFF function to parse GFF or GTF files and create a TxDb
txdb       <- makeTxDbFromGFF(file = file_annot, format = file_format)

# Use the readPeakFile function to read peaks files and generate an object
data_peaks <- readPeakFile(file_peaks)

# getPromoters retrieves the positions 3kb upstream and downstream of promoters
promoter   <- getPromoters(TxDb = txdb, upstream=3000, downstream=3000)

# getTagMatrix converts peaks results into a matrix with the window set to promoters
tagMatrix  <- getTagMatrix(data_peaks, windows = promoter)
peaksAnno  <- annotatePeak(data_peaks, TxDb=txdb, tssRegion=c(-3000, 3000), 
                           verbose=FALSE, addFlankGeneInfo=TRUE, flankDistance=5000)

# Convert to a table
annot_tab <- as.data.frame(peaksAnno)
colnames(annot_tab)[1] <- "#seqnames"

write.table(x = annot_tab, file = paste0(output_dir, "/", sample_name, ".peaks.gene_annotation.xls"),
            quote = FALSE, sep = "\t", row.names = FALSE)
#cat(paste0(paste(unique(annot_tab$geneId),collapse="\n"),"\n"), file=paste0(output_dir, "/", sample_name, ".peaks.gene_annotation.list"))

# Plot
# Visualize peak files: X-axis: chromosome positions; Y-axis: peaks positions and scores on each chromosome
pdf(file = paste0(output_dir, "/", sample_name, ".peaks.chromosome_distribution.pdf"))
covplot(data_peaks, weightCol = 'V5', xlab = "Chromosome Size (bp)", 
        title = paste0(sample_name, " Peaks over Chromosomes"))
dev.off()

# tagHeatmap function views the distribution of peaks in the 3kb region upstream and downstream of TSS (transcription start site) (heatmap)
pdf(file = paste0(output_dir, "/", sample_name, ".peaks.around_tss_heatmap.pdf"))
tagHeatmap(tagMatrix, xlab="Genomic Region (5'->3')", title = paste0(sample_name," peak distribution arround TSS"))
dev.off()

# plotAvgProf function plots a profile describing the average distribution
pdf(file = paste0(output_dir, "/", sample_name, ".peaks.around_tss_profile.pdf"), 6, 4)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", 
            ylab = "Peak Count Frequency", title = paste0(sample_name," peak distribution arround TSS"))
dev.off()

# plotDistToTSS function calculates the percentage of binding sites upstream and downstream of the nearest gene's TSS and visualizes their distribution
pdf(file = paste0(output_dir, "/", sample_name, ".peaks.around_tss_bar.pdf"), 7, 3)
plotDistToTSS(peaksAnno,title="Distribution of transcription factor-binding loci \n relative to TSS")
dev.off()

# plotAnnoBar function visualizes the composition of peaks: bar plot
pdf(file = paste0(output_dir, "/", sample_name, ".peaks.feature_distribution_bar.pdf"), 7, 3)
plotAnnoBar(x = peaksAnno, ylab = "Percentage(%)", title = paste0(sample_name, " Feature Distribution"))
dev.off()

# plotAnnoPie function visualizes the composition of peaks: pie chart
pdf(file = paste0(output_dir, "/", sample_name, ".peaks.feature_distribution_pie.pdf"), 6, 4)
plotAnnoPie(x = peaksAnno, pie3D = FALSE)
dev.off()

#__END__
