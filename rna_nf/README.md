# nf-core/atacseq

## 运行流程

```bash
# 本地运，如果能够访问github，可以直接运行：
nextflow run nf-core/atacseq -r 2.1.2 -name lw_atac_nf -profile singularity -resume -params-file nf-params.json
# 如果无法访问github，可以下载到本地运行。
# 首先在能够访问github的机器上下载nf-core/atacseq：
nf-core pipelines download atacseq -x tar.gz -s singularity
# 按照提示，设置好singularity缓存目录，例如：～/sifs。
# 然后等待下载完成，最后会得到一个tar.gz文件。
# 把这个文件上传到运行pipeline的机器上，同样需要设置好singularity缓存目录。
# 并且把下载pipeline时下载好的sigularity镜像缓存同样上传到运行pipeline的机器上的singularity缓存目录。
# 然后在运行pipeline的机器上的工作目录运行：
nextflow run nf-core-rnaseq_3.18.0/3_18_0/ -name lw_rna_nf -profile singularity -resume -params-file nf-params.json

# 合并bam文件
mamba activate samtools
samtools merge results/star_salmon/Ava.bam results/star_salmon/Ava_rep1.markdup.sorted.bam results/star_salmon/Ava_rep2.markdup.sorted.bam results/star_salmon/Ava_rep3.markdup.sorted.bam
samtools index results/star_salmon/Ava.bam
bamCoverage -b results/star_salmon/Ava.bam -o results/star_salmon/Ava.bw --extendReads --normalizeUsing RPKM

samtools merge results/star_salmon/DZ.bam results/star_salmon/DZ_rep1.markdup.sorted.bam results/star_salmon/DZ_rep2.markdup.sorted.bam results/star_salmon/DZ_rep3.markdup.sorted.bam
samtools index results/star_salmon/DZ.bam
bamCoverage -b results/star_salmon/DZ.bam -o results/star_salmon/DZ.bw --extendReads --normalizeUsing RPKM
```