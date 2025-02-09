# nf-core/atacseq

## Run Process

```bash
# Run it locally, if you have access to github, you can run it directly:
nextflow run nf-core/atacseq -r 2.1.2 -name lw_atac_nf -profile singularity -resume -params-file nf-params.json
# If you do not have access to github, you can download it locally.
# First, download the nf-core/atacseq pipeline to your local machine:
nf-core pipelines download atacseq -x tar.gz -s singularity
# Follow the instructions to set up the singularity cache directory, for example: ~/sifs.
# After the download is complete, you will have a tar.gz file.
# Upload this file to the machine where you will run the pipeline and set up the singularity cache directory there as well.
# Then, in the working directory on the machine where you will run the pipeline, run:
nextflow run nf-core-rnaseq_3.18.0/3_18_0/ -name lw_rna_nf -profile singularity -resume -params-file nf-params.json

# Merge bam files
mamba activate samtools
samtools merge results/star_salmon/Ava.bam results/star_salmon/Ava_rep1.markdup.sorted.bam results/star_salmon/Ava_rep2.markdup.sorted.bam results/star_salmon/Ava_rep3.markdup.sorted.bam
samtools index results/star_salmon/Ava.bam
bamCoverage -b results/star_salmon/Ava.bam -o results/star_salmon/Ava.bw --extendReads --normalizeUsing RPKM

samtools merge results/star_salmon/DZ.bam results/star_salmon/DZ_rep1.markdup.sorted.bam results/star_salmon/DZ_rep2.markdup.sorted.bam results/star_salmon/DZ_rep3.markdup.sorted.bam
samtools index results/star_salmon/DZ.bam
bamCoverage -b results/star_salmon/DZ.bam -o results/star_salmon/DZ.bw --extendReads --normalizeUsing RPKM
```