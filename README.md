```
conda create -n sgc-arg snakemake-minimal mamba
conda activate sgc-arg

# on farm:
snakemake -j 4 --use-conda --rerun-incomplete --restart-times 1 --latency-wait 15 --resources mem_mb=400000 --default-resources runtime=1440 --cluster "sbatch -t {resources.runtime} -J ibd_smake -p medium -n 1 -N 1 -c {threads} --mem={resources.mem_mb}"
```
