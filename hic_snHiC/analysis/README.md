```bash
# loops gene analysis
cat 09_Loop_detection_HiCexplorer/DZ/DZ_mapQ15_10kb_loops.bedpe | awk '{print $1"\t"$2"\t"$3}' > analysis/loops_gene.bed
cat 09_Loop_detection_HiCexplorer/DZ/DZ_mapQ15_10kb_loops.bedpe | awk '{print $4"\t"$5"\t"$6}' >> analysis/loops_gene.bed
cat analysis/loops_gene.bed | sort -k1,1 -k2,2n | uniq > analysis/loops_gene_sorted.bed
# activate conda environment
mamba activate bedtools
bedtools intersect -a analysis/loops_gene_sorted.bed \
-b /home/gyk/reference/gff/Bomo_genes_chr.bed \
-wa -wb | awk '{print $8}' | sort | uniq > analysis/loops_geneid.list

# activate conda environment
mamba activate R
Rscript analysis/enrichx.R -i analysis/loops_gene
id.list -o analysis/results/enrichment_results -d org.Bmori.eg.db
```

## Differential Loops Analysis

```bash
awk 'FNR==NR {key=$1 FS $2 FS $3 FS $4 FS $5 FS $6; seen[key]++; next} {key=$1 FS $2 FS $3 FS $4 FS $5 FS $6; if (seen[key] != 1) print}' 09_Loop_detection_HiCexplorer/Ava/Ava_mapQ15_10kb_loops.bedpe 09_Loop_detection_HiCexplorer/DZ/DZ_mapQ15_10kb_loops.bedpe > 09_Loop_detection_HiCexplorer/DZ_Ava_10kb_loops.bedpe

cat 09_Loop_detection_HiCexplorer/DZ_Ava_10kb_loops.bedpe | awk '{print $1"\t"$2"\t"$3}' > analysis/diff_loops_gene.bed
cat 09_Loop_detection_HiCexplorer/DZ_Ava_10kb_loops.bedpe | awk '{print $4"\t"$5"\t"$6}' >> analysis/diff_loops_gene.bed

mamba activate bedtools
bedtools intersect -a analysis/diff_loops_gene.bed \
-b /home/gyk/reference/gff/Bomo_genes_chr.bed \
-wa -wb | awk '{print $8}' | sort | uniq > analysis/diff_loops_geneid.list

# annotation
grep -Ff analysis/diff_loops_geneid.list analysis/Bmo_annotations.tsv | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4}' > analysis/diff_loops_gene_annotation.tsv

```

## Plot Differential Loops

```bash
mamba activate hicexplorer
hicPlotMatrix -m 06_Interaction_matrices_normalized_and_corrected/corrected_matrices/Ava/cool_format/Ava_mapQ15_10kb_normalized_corrected.cool \
-o analysis/ava_chr12_plot2.pdf --dpi 300 --log1p \
--region chr12:10000000-12000000 \
--loops 09_Loop_detection_HiCexplorer/Ava/Ava_mapQ15_10kb_loops.bedpe

hicPlotMatrix -m 06_Interaction_matrices_normalized_and_corrected/corrected_matrices/DZ/cool_format/DZ_mapQ15_10kb_normalized_corrected.cool \
-o analysis/dz_chr12_plot2.pdf --dpi 300 --log1p \
--region chr12:10000000-12000000 \
--loops 09_Loop_detection_HiCexplorer/DZ/DZ_mapQ15_10kb_loops.bedpe

hicPlotMatrix -m 06_Interaction_matrices_normalized_and_corrected/corrected_matrices/Ava/cool_format/Ava_mapQ15_10kb_normalized_corrected.cool \
-o analysis/chr18_plot.pdf --dpi 300 --log1p \
--region chr18:9000000-13000000 \
--loops 09_Loop_detection_HiCexplorer/Ava/Ava_mapQ15_10kb_loops.bedpe
```

```bash
# Differential Contacts Gene Analysis
awk 'NR > 1 && $8 > 0' 12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/Ava_vs_DZ/Ava_vs_DZ_10kb_SELFISH.txt > analysis/selfish_filtered_result.txt

cat analysis/selfish_filtered_result.txt | awk '{print $1"\t"$2"\t"$3}' > analysis/diff_contacts_gene.bed
cat analysis/selfish_filtered_result.txt | awk '{print $4"\t"$5"\t"$6}' >> analysis/diff_contacts_gene.bed
cat analysis/diff_contacts_gene.bed | sort -k1,1 -k2,2n | uniq | awk '{print "chr"$1"\t"$2"\t"$3}' > analysis/diff_contacts_gene_sorted.bed

mamba activate bedtools
bedtools intersect -a analysis/diff_contacts_gene_sorted.bed \
-b /home/gyk/project/lw_atac_nf/analysis/results/DZ_Ava_diff_consensus_peaks.bed -F 0.5 \
-wa -wb > analysis/diff_contacts_peak.bed

bedtools intersect -a analysis/diff_contacts_gene_sorted.bed \
-b /home/gyk/project/lw_atac_nf/analysis/results/DZ_Ava_diff_consensus_peaks.bed -F 0.5 \
-wa -wb | awk '{print $8}' | sort | uniq > analysis/diff_contacts_peak_id.list

# annotation
grep -Ff analysis/diff_contacts_geneid.list analysis/Bmo_annotations.tsv | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4}' > analysis/diff_contacts_gene_annotation.tsv

```

## Find Genes with Significant Changes in Loops and Expression Levels

```bash
# Extract the first six columns and sort them
awk '{print $1, $2, $3, $4, $5, $6}' 09_Loop_detection_HiCexplorer/Ava/Ava_mapQ15_10kb_loops.bedpe | sort > analysis/ava_loops.txt
awk '{print $1, $2, $3, $4, $5, $6}' 09_Loop_detection_HiCexplorer/DZ/DZ_mapQ15_10kb_loops.bedpe | sort > analysis/dz_loops.txt

# Compare files using comm
# Only in file1
comm -23 analysis/ava_loops.txt analysis/dz_loops.txt > analysis/only_in_ava.txt

# Only in file2
comm -13 analysis/ava_loops.txt analysis/dz_loops.txt > analysis/only_in_dz.txt

# only in ava genes
cat analysis/only_in_ava.txt | awk '{print $1"\t"$2"\t"$3}' > analysis/only_in_ava_loops_gene.bed
cat analysis/only_in_ava.txt | awk '{print $4"\t"$5"\t"$6}' >> analysis/only_in_ava_loops_gene.bed

mamba activate bedtools
bedtools intersect -a analysis/only_in_ava_loops_gene.bed \
-b /home/gyk/reference/gff/Bomo_genes_chr.bed \
-F 0.5 -wa -wb | awk '{print $8}' | sort | uniq > analysis/only_in_ava_loops_geneid.list

# only in dz genes
cat analysis/only_in_dz.txt | awk '{print $1"\t"$2"\t"$3}' > analysis/only_in_dz_loops_gene.bed
cat analysis/only_in_dz.txt | awk '{print $4"\t"$5"\t"$6}' >> analysis/only_in_dz_loops_gene.bed

mamba activate bedtools
bedtools intersect -a analysis/only_in_dz_loops_gene.bed \
-b /home/gyk/reference/gff/Bomo_genes_chr.bed \
-F 0.5 -wa -wb | awk '{print $8}' | sort | uniq > analysis/only_in_dz_loops_geneid.list

grep -Ff analysis/only_in_ava_loops_geneid_up.list analysis/Bmo_annotations.tsv | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4}' > analysis/only_in_ava_loops_geneid_up_annotation.tsv

# Create a file for differential links
awk 'NR > 1 && $8 > 0 {print "chr"$1"\t"$2"\t"$3"\tchr"$4"\t"$5"\t"$6"\t"$8}' 12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/Ava_vs_DZ/Ava_vs_DZ_10kb_SELFISH.txt > analysis/Ava_link.tsv
awk 'NR > 1 && $8 < 0 {print "chr"$1"\t"$2"\t"$3"\tchr"$4"\t"$5"\t"$6"\t"($8 < 0 ? -$8 : $8)}' 12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/Ava_vs_DZ/Ava_vs_DZ_10kb_SELFISH.txt > analysis/DZ_link.tsv