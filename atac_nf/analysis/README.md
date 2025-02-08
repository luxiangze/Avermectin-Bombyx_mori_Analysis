# Analysis script
## Filtering consensus peaks

```bash
# activate conda environment
mamba activate pandas
python3 analysis/filter_consensus_peaks.py results/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.boolean.txt results/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.annotatePeaks.txt > analysis/consensus_peaks.mLb.clN.filtered.txt
```

## plot merged replicate distance of peaks to TSS

Method 1:Plot one or more samples

```bash
# activate conda environment and make directory
mkdir -p analysis/merged_replicate
mamba activate deeptools
# create matrix file
computeMatrix reference-point --referencePoint TSS \
    -R results/genome/Bomo_gene_models_chr.bed \
    -S results/bwa/merged_replicate/bigwig/DZ.mRp.clN.bigWig results/bwa/merged_replicate/bigwig/Ava.mRp.clN.bigWig \
    --samplesLabel DZ Ava \
    -p 40 -b 3000 -a 3000 --skipZeros \
    --outFileName analysis/results/merged_replicate/DZ-Ava.matrix_TSS.cpm.gz \
    --outFileNameMatrix analysis/results/merged_replicate/DZ-Ava.matrix_TSS.cpm.tab \
    --outFileSortedRegions analysis/results/merged_replicate/DZ-Ava.matrix_TSS.cpm.bed

# plot heatmap
plotHeatmap --matrixFile analysis/results/merged_replicate/DZ-Ava.matrix_TSS.cpm.gz \
    --outFileName analysis/results/merged_replicate/DZ-Ava.heatmap_TSS.cpm.pdf \
    --colorMap RdBu \
	--heatmapWidth 6 --heatmapHeight 15 \
	--whatToShow 'plot, heatmap and colorbar' \
	--legendLocation best --refPointLabel TSS \
	--xAxisLabel 'Distance to TSS (bp)' --yAxisLabel 'CPM' \
	--regionsLabel '' --plotFileFormat pdf --dpi 720

# plot profile
plotProfile --matrixFile analysis/results/merged_replicate/DZ-Ava.matrix_TSS.cpm.gz \
	--outFileName analysis/results/merged_replicate/DZ-Ava.profile_TSS.cpm.pdf \
	--refPointLabel TSS --yAxisLabel 'CPM' --perGroup \
	--regionsLabel 'treat and input signal around TSS' \
	--plotFileFormat pdf --dpi 720
```

# Method 2:Compare with the input (control) sample (enrichment in the gene body region).

```bash
bamCompare -b1 results/bwa/merged_replicate/Ava.mRp.clN.sorted.bam -b2 results/bwa/merged_replicate/DZ.mRp.clN.sorted.bam \
	--operation log2 --normalizeUsing CPM --scaleFactorsMethod None --binSize 50 \
	-o analysis/results/merged_replicate/DZ-Ava.cpm_log2ratio.bw --numberOfProcessors 60

computeMatrix scale-regions -R results/genome/Bomo_gene_models_chr.bed \
	-S analysis/results/merged_replicate/DZ-Ava.cpm_log2ratio.bw --samplesLabel DZ-Ava -p 60 -b 3000 -a 3000 \
	--skipZeros --outFileName analysis/results/merged_replicate/DZ-Ava.matrix_gene.cpm_log2ratio.gz \
	--outFileNameMatrix analysis/results/merged_replicate/DZ-Ava.matrix_gene.cpm_log2ratio.tab \
	--outFileSortedRegions analysis/results/merged_replicate/DZ-Ava.matrix_gene.cpm_log2ratio.bed

plotHeatmap --matrixFile analysis/results/merged_replicate/DZ-Ava.matrix_gene.cpm_log2ratio.gz \
	--outFileName analysis/results/merged_replicate/DZ-Ava.heatmap_gene.cpm_log2ratio.pdf \
	--heatmapWidth 6 --heatmapHeight 15 \
	--colorMap RdBu --whatToShow 'plot, heatmap and colorbar' \
	--legendLocation best --refPointLabel gene --xAxisLabel 'Distance to gene(bp)' \
	--yAxisLabel 'log2(Ava/DZ)' --regionsLabel '' --plotFileFormat pdf --dpi 720

plotProfile --matrixFile analysis/results/merged_replicate/DZ-Ava.matrix_gene.cpm_log2ratio.gz \
	--outFileName analysis/results/merged_replicate/DZ-Ava.profile_gene.cpm_log2ratio.pdf \
	--refPointLabel gene --yAxisLabel 'log2(Ava/DZ)' --perGroup \
	--regionsLabel 'DZ-Ava signal around gene' --plotFileFormat pdf --dpi 720
```

## diff peak analysis

```bash
cat analysis/results/DZ_Ava_diff_consensus_peaks.bed | awk '{print $5}' |sort |uniq > analysis/results/diff_peak_gene_id.txt

# enrichment
mamba activate R
Rscript analysis/enrichx.R -i analysis/results/diff_peak_gene_id.txt -o analysis/results/enrichment_results -d org.Bmori.eg.db

```

