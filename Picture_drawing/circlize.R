# 设置工作目录
setwd("/home/gyk/project/lw")

start_time <- Sys.time()

# 导入circlize、ggplot2、grid和ggplotify包
library(circlize)
library(ggplot2)
library(grid)
library(ggplotify)

pdf("circlize_plot.pdf", width = 8, height = 8)

# legend

p1 <- ggplot() +
  # Inner and outer labels
  geom_text(aes(x = 1, y = c(2, 9),
                label = c("Inner", "Outer")),
            size = 3) +
  # Main vertical arrow
  geom_segment(aes(x = 2, y = 10, xend = 2, yend = 0),
               arrow = arrow(angle = 20, type = "closed", length = unit(5, "mm")),
               linewidth = 1) +
  # # Vertical line for categories
  # geom_segment(aes(x = 2, y = 10, xend = 2, yend = 1), linewidth = 1) +
  # color points
  geom_point(aes(x = 3, y = c(4, 6, 8)),
             shape = 15, color = c("red", "purple", "blue"), size = 3) +
  # Main text labels
  geom_text(aes(x = c(2.5, 3.25, 3.25, 3.25), y = c(2, 4, 6, 8), hjust = 0,
                label = c("Genome Compartment", "ATAC", "RNA", "Gene")),
            size = 3) +
  # Vertical segments
  geom_segment(aes(x = c(2.5, 4.5), y = c(1, 1),
                   xend = c(2.5, 4.5), yend = c(1.5, 1.5)),
               linewidth = 0.5) +
  # Horizontal segments
  geom_segment(aes(x = 2.5, y = 1.5,
                   xend = 4.5, yend = 1.5),
               linewidth = 0.5) +
  # color points
  geom_point(aes(x = c(2.5, 4.5), y = 1),
             shape = 15, color = c("salmon", "skyblue"), size = 3) +
  # Segment values
  geom_text(aes(x = c(2.75, 4.75), y = 1,
                label = c("A", "B")),
            size = 3, hjust = 0) +
  # Remove background and axes
  theme_void()
# Display the plot
# p1


# 创建染色体信息
chrom_df <- read.table("/home/gyk/reference/genome/Bomo_genome_assembly.chrom.bed", sep = "\t", colClasses = c("character", "numeric", "numeric"))
colnames(chrom_df) <- c("chr", "start", "end")

# 初始化环形图参数
print("Initializing circos plot")

circos.par(
  "track.height" = 0.05,
  "start.degree" = 90, # 调整起始角度为90度
  "gap.degree" = 2
) # 设置染色体间间隔为2度

# 创建环形图基础层
circos.initializeWithIdeogram(chrom_df,
  plotType = c("axis", "labels")
)

print("Adding gene density tracks")

# 读取并处理基因密度数据
gene_density <- read.table("/home/gyk/reference/genome/Bomo_gene_density.txt", header = TRUE)

# 绘制基因密度轨道
circos.genomicTrack(
  gene_density,
  track.height = 0.05,
  bg.border = NA,
  bg.col = "#FFFFFF",
  ylim = c(0, 1), # 明确设置y轴范围
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value,
      ytop.column = 1,
      ybottom = 0,
      col = ifelse(value[[2]] == "gene", "blue", "white"), # 深蓝色表示基因区域
      border = NA
    )
  }
)

print("Adding RNA-seq tracks")

# 添加RNA-seq数据
gene_bed <- read.table("/home/gyk/project/lw_rna_nf/results/star_salmon/stringtie/DZ.gene.abundance.txt", header = TRUE, sep = "\t")
gene_bed  <- gene_bed[, c("Reference", "Start", "End", "TPM")]

# 计算RNA-seq数据的值范围
summary(gene_bed$TPM)

# 进行log2转换
gene_bed$TPM <- log2(gene_bed$TPM + 1)
rna_max <- max(gene_bed$TPM)
summary(gene_bed$TPM)

circos.genomicTrack(gene_bed,
  track.height = 0.1,
  bg.border = NA,
  bg.col = "#FFFFFF",
  ylim = c(0, rna_max), # 使用实际的数据范围
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value,
      type = "h",
      col = "purple",
      lwd = 0.5
    )
  }
)

print("Adding ATAC-seq tracks")

# 添加ATAC-seq数据
atac_df <- read.table("/home/gyk/project/lw_atac_nf/results/bwa/merged_replicate/macs2/narrow_peak/DZ.mRp.clN_peaks.xls", comment.char = "#", header = TRUE)
atac_df <- atac_df[, c(1, 2, 3)]
colnames(atac_df) <- c("chr", "start", "end")

# 设置分辨率为 100kb
resolution <- 100000

# 加载所需的包
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")

library(dplyr)
library(tidyr)

# 创建一个函数来统计每条染色体的峰数量
count_peaks <- function(data, resolution) {
  data %>%
    mutate(bin_start = floor(start / resolution) * resolution,
           bin_end = floor(end / resolution) * resolution) %>%
    # 为了考虑跨区间的情况，将每个峰展开成所有跨越的 bin
    rowwise() %>%
    mutate(bins = list(seq(bin_start, bin_end, by = resolution))) %>%
    unnest(cols = c(bins)) %>%
    group_by(chr, bins) %>%
    summarise(peak_count = n(), .groups = "drop")
}

# 统计峰数量
peak_counts <- count_peaks(atac_df, resolution)

# 查看结果
print(peak_counts)

# 统计ATAC-seq数据的值范围
summary(atac_df["log10(qvalue)"])
atac_df["log10(qvalue)"] <- log2(atac_df["log10(qvalue)"] + 1)

atac_max <- max(atac_df["log10(qvalue)"])
summary(atac_df["log10(qvalue)"])

circos.genomicTrack(atac_df,
  track.height = 0.1,
  bg.border = NA,
  bg.col = "#FFFFFF",
  ylim = c(0, atac_max), # 使用实际的数据范围
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value,
      type = "h",
      col = "#FF0000",
      lwd = 0.5
    )
  }
)

print("Adding compartment tracks")

# 添加基因组compartment数据
compartment <- read.table("/home/gyk/project/lw_hic_snHiC/10_Compartments_detection_dcHiC/20kb_resolution/DifferentialResult/all_vs_all/viz/files_compartment_beds/intra_DZ_mapQ15_20kb_PC_compartments_sorted.bed", skip = 1, header = FALSE)
compartment <- compartment[, c("V1", "V2", "V3", "V4", "V5")]
colnames(compartment) <- c("chr", "start", "end", "compartment", "score")

# 绘制compartment轨道,A区室绘制成浅红色,B区绘制成浅蓝色
circos.genomicTrack(compartment,
  track.height = 0.05,
  bg.border = NA,
  bg.col = "#FFFFFF",
  ylim = c(0, 1), # 明确设置y轴范围
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value,
      ytop = 1,
      ybottom = 0,
      col = ifelse(value$compartment == "A", "salmon", "skyblue"),
      border = NA
    )
  }
)

# adding legend
p2 <- as.grob(p1)
vp <- viewport(x = .5, y = .5, width = .3, height = .3)
pushViewport(vp)
grid.draw(p2)
upViewport()

# 清除环形图设置
circos.clear()
dev.off()

end_time <- Sys.time()
time <- end_time - start_time
print(time)

print("Done")