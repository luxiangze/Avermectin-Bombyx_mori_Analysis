# Dry-run and run

```bash
# activate conda environment
mamba activate snHiC

# Dry-run
snakemake \
-s /home/gyk/project/lw_hic_snHiC/snHiC.snakefile \
--configfile /home/gyk/project/lw_hic_snHiC/snHiC_config.yaml \
--cores 30 \
-n

# Run
snakemake \
-s /home/gyk/project/lw_hic_snHiC/snHiC.snakefile \
--configfile /home/gyk/project/lw_hic_snHiC/snHiC_config.yaml \
--cores 30
```