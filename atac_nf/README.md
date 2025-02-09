# nf-core/atacseq

## Workflow

```bash
# Local execution: If you can access GitHub, you can directly run:
nextflow run nf-core/atacseq -r 2.1.2 -name lw_atac_nf -profile singularity -resume -params-file nf-params.json
# If you cannot access GitHub, you can download and run it locally.
# First, download nf-core/atacseq on a machine that can access GitHub:
nf-core pipelines download atacseq -x tar.gz -s singularity
# Follow the prompts to set up the Singularity cache directory, e.g., ~/sifs.
# Then wait for the download to complete, and you will get a tar.gz file.
# Upload this file to the machine running the pipeline, and also set up the Singularity cache directory.
# Also upload the Singularity image cache downloaded during the pipeline download to the Singularity cache directory on the machine running the pipeline.
# Then, in the working directory on the machine running the pipeline, run:
nextflow run nf-core-atacseq_2.1.2/2_1_2/  -name lw_atac_nf -profile singularity -resume -params-file nf-params.json
```