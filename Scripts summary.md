### Submit jobs
```bash
# main
snakemake -s snakefile_concoct1.py -j 1 --cluster "sbatch -A naiss2024-22-820 -p main -N 1 -c 128 -t 11:00:00 -J concoct --mail-user ssamantha.data@gmail.com --mail-type=ALL" --keep-incomplete

# shared
snakemake -s abundance.single.change.smk -j 1 --cluster "sbatch -A naiss2024-22-1134 -p shared --ntasks=1 --cpus-per-task=20 -t 02:00:00 -J ecGEM --mail-user ssamantha.data@gmail.com --mail-type=ALL" --keep-incomplete
```

### Snakemake-MetaGEM
1. Trimming raw reads (fastp): `snakefile_megahit.smk`
2. Assembling raw reads (megahit): `megahit.smk`
3. Refining contigs (Metawrap) \
   3.1 Crossmap: `crossmap.smk` \
   3.2.1 Maxbin: `snakefile_maxbin1.smk` need to double confirm the parameter is 0.8 or 0.9. \
   3.2.2 Metabat2: \
   3.2.3 Concoct: `concoct.smk` \
   3.3 Refine: `refine.smk` \
   3.4 Reassemble: `reassemble.smk`
4. Extracting protein and DNA bins: `extract.smk`
5. Profiling taxonomy (GTDB-tk): `taxo.smk`
6. Creating GEMs (CarveMe): `carveme.smk`
7. Converting GEMs to ecGEMs:`GEMtoECGEM.smk`
8. Calculating abundance (bwa):  `abundance.smk`
9. Annotating prokaryotic genome (Prokka): ``

### MatLab-scripts
1. Extracting metabolites' name from GEMs.
2. Converting GEMs to ecGEMs and FVA.
3. Changed functions of GECKO3.
4. Other related files: `metName.inhouseDB.txt` 

### R-scripts 
1. Extracting species KEGG ID and Uniprot ID:
2. Extracting SMILEs (PubChem API)
3. Metagenomics analysis
4. Constructing Metabolite-Microbe network
5. Reporter metabolite analysis
