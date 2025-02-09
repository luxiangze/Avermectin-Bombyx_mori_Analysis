# nf-core/differentialabundance

## Workflow

```bash
# Local execution: If you can access GitHub, you can directly run:
nextflow run nf-core/differentialabundance -r 1.5.0 -name lw_DGE_nf -profile singularity -resume -params-file nf-params.json
# If you cannot access GitHub, you can download and run it locally.
# First, download nf-core/differentialabundance on a machine that can access GitHub:
nf-core pipelines download differentialabundance -x tar.gz -s singularity
# Follow the prompts to set up the Singularity cache directory, e.g., ~/sifs.
# Then wait for the download to complete, and you will get a tar.gz file.
# Upload this file to the machine running the pipeline, and also set up the Singularity cache directory.
# Also upload the Singularity image cache downloaded during the pipeline download to the Singularity cache directory on the machine running the pipeline.
# Then, in the working directory on the machine running the pipeline, run:
nextflow run nf-core-differentialabundance_1.5.0/1_5_0 -name lw_rnaseq_dge -profile docker -resume -params-file nf-params.json
```

# plot
`analysis/DGE_analysis.r` was used to plot volcano and heatmap.`enrichx.R` was used to plot GO and KEGG enrichment.

```bash
mamba activate R
Rscript analysis/enrichx.R -i analysis/results/DEG_id.list -o analysis/results/enrichment_results/ -d org.Bmori.eg.db
```