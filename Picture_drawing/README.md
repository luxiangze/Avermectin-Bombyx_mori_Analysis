## make tracks file

```bash
# Activate the pyGenomeTracks environment
mamba activate pgt

# Generate tracks file
make_tracks_file \
--trackFiles \
"/home/gyk/project/lw_hic_snHiC/06_Interaction_matrices_normalized_and_corrected/corrected_matrices/Ava/h5_format/Ava_mapQ15_20kb_normalized_corrected.h5" \
"/home/gyk/project/lw_hic_snHiC/07_TADs_calling_HiCexplorer/Ava/20kb_resolution/Ava_mapQ15_20kb_domains.bed" \
"/home/gyk/project/lw_hic_snHiC/06_Interaction_matrices_normalized_and_corrected/corrected_matrices/DZ/h5_format/DZ_mapQ15_20kb_normalized_corrected.h5" \
"/home/gyk/project/lw_hic_snHiC/07_TADs_calling_HiCexplorer/DZ/20kb_resolution/DZ_mapQ15_20kb_domains.bed" \
"/home/gyk/project/lw_rna_nf/results/star_salmon/bigwig/Ava.forward.bw" \
"/home/gyk/project/lw_rna_nf/results/star_salmon/bigwig/DZ.forward.bw" \
"/home/gyk/project/lw_atac_nf/results/bwa/merged_replicate/macs2/narrow_peak/Ava.mRp.clN_peaks.narrowPeak" \
"/home/gyk/project/lw_atac_nf/results/bwa/merged_replicate/macs2/narrow_peak/DZ.mRp.clN_peaks.narrowPeak" \
"/home/gyk/reference/gtf/Bomo_gene_models_chr.gtf" \
-o templet_tracks.ini
```

## plot tracks

```bash
# Generate tracks file
pyGenomeTracks --tracks tracks.ini --region chr12:6167018-6382366 -o KWMTBOMO07140_image.pdf --dpi 300

for i in `cat gene.list`
do
    echo $i
    region=`grep $i /home/gyk/reference/gtf/Bomo_gene_models_chr.gtf | grep -P '\tgene\t' | awk '{print $1":"$4 - 20000"-"$5 + 20000}'`
    # echo $region
    pyGenomeTracks --tracks tracks.ini --region $region -o ${i}_tracks.png --dpi 300 -t $i
done

# hic
pyGenomeTracks --tracks hic_tracks.ini --region chr12:4767018-7782366 -o KWMTBOMO07140_hic_image.pdf --dpi 300

# loops
hicPlotMatrix -m /home/gyk/project/lw_hic_snHiC/06_Interaction_matrices_normalized_and_corrected/corrected_matrices/Ava/h5_format/Ava_mapQ15_20kb_normalized_corrected.h5 -o plot.png --region chr13:4376570-4607580 --loops /home/gyk/project/lw_hic_snHiC/09_Loop_detection_HiCexplorer/Ava/Ava_mapQ15_20kb_loops.bedpe
```

## gene annotation

```bash
grep -Ff gene.list /home/gyk/project/lw_hic_snHiC/analysis/Bmo_annotations.tsv | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4}' > gene.list.annotation.tsv
```

## Batch generation track

```bash
for i in `cat gene.list`
do
    echo $i
    hic_region=`grep $i /home/gyk/project/lw_hic_snHiC/analysis/Bomo_genes_chr.bed | awk '{print $1":"$2 - 1500000"-"$3 + 1500000}'`
    pyGenomeTracks --tracks hic_tracks.ini --region $hic_region -o hic/${i}_hic_tracks.pdf --dpi 300
    region=`grep $i /home/gyk/project/lw_hic_snHiC/analysis/Bomo_genes_chr.bed | awk '{print $1":"$2 - 100000"-"$3 + 100000}'`
    # echo $region
    pyGenomeTracks --tracks tracks.ini --region $region -o hic/${i}_tracks.pdf --dpi 300
done

# links 
# If the 8 column is greater than 0, it is written to link_increase.txt, if it is less than 0, it is written to link_decrease.txt
awk '{if($8>0){print "chr"$1"\t"$2"\t"$3"\tchr"$4"\t"$5"\t"$6"\t"$7"\t"$8}}' /home/gyk/project/lw_hic_snHiC/12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/Ava_vs_DZ/Ava_vs_DZ_1kb_SELFISH.txt > link_increase.txt
awk '{if($8<=0){print "chr"$1"\t"$2"\t"$3"\tchr"$4"\t"$5"\t"$6"\t"$7"\t"$8}}' /home/gyk/project/lw_hic_snHiC/12_Grouped_analyses/G_Differential_contacts_analyses_SELFISH/Ava_vs_DZ/Ava_vs_DZ_1kb_SELFISH.txt > link_decrease.txt

```