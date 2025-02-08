# nf-core/differentialabundance

## 运行流程

```bash
# 本地运，如果能够访问github，可以直接运行：
nextflow run nf-core/differentialabundance -r 1.5.0 -name lw_DGE_nf -profile singularity -resume -params-file nf-params.json
# 如果无法访问github，可以下载到本地运行。
# 首先在能够访问github的机器上下载nf-core/differentialabundance：
nf-core pipelines download differentialabundance -x tar.gz -s singularity
# 按照提示，设置好singularity缓存目录，例如：～/sifs。
# 然后等待下载完成，最后会得到一个tar.gz文件。
# 把这个文件上传到运行pipeline的机器上，同样需要设置好singularity缓存目录。
# 并且把下载pipeline时下载好的sigularity镜像缓存同样上传到运行pipeline的机器上的singularity缓存目录。
# 然后在运行pipeline的机器上的工作目录运行：
nextflow run nf-core-differentialabundance_1.5.0/1_5_0 -name lw_rnaseq_dge -profile docker -resume -params-file nf-params.json
```

# plot
`analysis/DGE_analysis.r` was used to plot volcano and heatmap.`enrichx.R` was used to plot GO and KEGG enrichment.

```bash
mamba activate R
Rscript analysis/enrichx.R -i analysis/results/DEG_id.list -o analysis/results/enrichment_results/ -d org.Bmori.eg.db
```