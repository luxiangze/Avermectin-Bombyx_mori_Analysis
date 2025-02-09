# Set working directory
setwd("/home/gyk/project/lw")

start_time <- Sys.time()

# Import circlize, ggplot2, grid, and ggplotify packages
library(circlize)
library(ggplot2)
library(grid)
library(ggplotify)

pdf("circlize_hic_plot.pdf", width = 8, height = 8)

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
  geom_point(aes(x = 2.75, y = c(4, 6, 8)),
             shape = 15, color = c("red", "black", "blue"), size = 3) +
  # Main text labels
  geom_text(aes(x = c(2.5, 3, 3, 3), y = c(2, 4, 6, 8), hjust = 0,
                label = c("Genome Compartment", "Hi-C Loop", "TAD Boundary", "Gene Density")),
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


# Create chromosome information
chrom_df <- read.table("/home/gyk/reference/genome/Bomo_genome_assembly.chrom.bed", sep = "\t", colClasses = c("character", "numeric", "numeric"))
colnames(chrom_df) <- c("chr", "start", "end")

# Initialize circular plot parameters
print("Initializing circos plot")

circos.par(
  "track.height" = 0.05,
  "start.degree" = 90, # Set start angle to 90 degrees
  "gap.degree" = 2
) # Set chromosome spacing to 2 degrees

# Create base layer for circular plot
circos.initializeWithIdeogram(chrom_df,
  plotType = c("axis", "labels")
)

print("Adding gene density tracks")

# Read and process gene density data
gene_density <- read.table("/home/gyk/reference/genome/Bomo_gene_density_100kb.bed", header = TRUE)

# Plot gene density track
circos.genomicTrack(
  gene_density,
  track.height = 0.1,
  bg.border = NA,
  bg.col = "#FFFFFF",
  ylim = c(0, 1), # Explicitly set y-axis range
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value,
      type = "h",
      col = "blue",
      lwd = 0.5
    )
  }
)

print("Adding TAD tracks")

# Add TAD data
tad_bed <- read.table("/home/gyk/project/lw_hic_snHiC/07_TADs_calling_HiCexplorer/DZ/100kb_resolution/DZ_mapQ15_100kb_boundaries.bed", header = FALSE, sep = "\t")
tad_bed  <- tad_bed[, c("V1", "V2", "V3")]

circos.genomicTrack(tad_bed,
  track.height = 0.1,
  bg.border = NA,
  bg.col = "#FFFFFF",
  ylim = c(0, 1),
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value,
      ybottom = 0,
      ytop = 1,
      type = "h",
      col = "black",
      lwd = 0.05
    )
  }
)

print("Adding loop link tracks")

# Add loop data
loop <- read.table("/home/gyk/project/lw_hic_snHiC/09_Loop_detection_HiCexplorer/DZ/DZ_mapQ15_20kb_loops.bedpe", header = FALSE, sep = "\t")
colnames(loop) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "pvalue")
loops1 <- loop[, c("chr1", "start1", "end1", "pvalue")]
colnames(loops1) <- c("chr", "start", "end", "pvalue")
loops2 <- loop[, c("chr2", "start2", "end2", "pvalue")]
colnames(loops2) <- c("chr", "start", "end", "pvalue")
# Merge two data frames
loops <- rbind(loops1, loops2)

# Calculate loop p-value data
loops$pvalue <- -log10(loops$pvalue)
summary(loops$pvalue)
p_max <- max(loops$pvalue)

circos.genomicTrack(loops,
  track.height = 0.1,
  bg.border = NA,
  bg.col = "#FFFFFF",
  ylim = c(0, p_max), # Use actual data range
  panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value,
      type = "h",
      col = "#FF0000",
      lwd = 0.5
    )
  }
)

print("Adding compartment tracks")

# Add genome compartment data
compartment <- read.table("/home/gyk/project/lw_hic_snHiC/10_Compartments_detection_dcHiC/20kb_resolution/DifferentialResult/all_vs_all/viz/files_compartment_beds/intra_DZ_mapQ15_20kb_PC_compartments_sorted.bed", skip = 1, header = FALSE)
compartment <- compartment[, c("V1", "V2", "V3", "V4", "V5")]
colnames(compartment) <- c("chr", "start", "end", "compartment", "score")

# Plot compartment track, A compartment in light red, B compartment in light blue
circos.genomicTrack(compartment,
  track.height = 0.05,
  bg.border = NA,
  bg.col = "#FFFFFF",
  ylim = c(0, 1), # Explicitly set y-axis range
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

# Clear circular plot settings
circos.clear()
dev.off()

end_time <- Sys.time()
time <- end_time - start_time
print(time)

print("Done")