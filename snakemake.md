1. Submit jobs
```bash
snakemake -s snakefile_concoct1.py -j 1 --cluster "sbatch -A naiss2024-22-820 -p main -N 1 -c 128 -t 11:00:00 -J concoct --mail-user ssamantha.data@gmail.com --mail-type=ALL" --keep-incomplete
```
