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
   3.2.1 Maxbin: `maxbin.smk` 0.9. \
   3.2.2 Metabat2: `metabat.smk` \
   3.2.3 Concoct: `concoct.smk` \
   3.3 Refine: `refine.smk` \
   3.4 Reassemble: `reassemble.smk`
4. Extracting protein and DNA bins: `extract.smk`
5. Profiling taxonomy (GTDB-tk): `taxo.smk`
6. Creating GEMs (CarveMe): `carveme.smk`
7. Converting GEMs to ecGEMs:`GEMtoECGEM.smk`
8. Calculating abundance (bwa):  `abundance.smk`
9. Annotating prokaryotic genome (Prokka): `prokka.smk`

### MatLab-scripts
1. Extracting metabolites' name from GEMs: `meta-gecko-extractMet.m`
2. Converting GEMs to ecGEMs and FVA: `meta-gecko.m`
3. Adapter template for customizing adapter files for each bin: `adapterTemplate.m`
4. Changed functions of GECKO3:
   1) `makeEcModel.m`: key function to change GEM to ecGEM.
   2) `getECfromDatabase_uniprot.m`
   3) `getECfromDatabase_kegg.m`
   4) `fuzzyKcatMatching.m`
   5) `loadConventionalGEM.m`: importModel set to false to allow exchange reations showing. 
6. Other related files: 
   1) `metName.inhouseDB.txt`: Extracted SMILEs by metabolites' names from PubChem API.
   2) `organism.txt`: KEGG prokaryote organism names and ID, for extracting the KEGG ID. 
   3) `taxonomy.R`: Fuzzy matching taxonomies for each bin to KEGG and Uniprot database to extract the KEGG ID and Uniprot ID. 
   4) `UniprotID1.txt` and `UniprotID2.txt`: Uniprot ID for extracting EC numbers from Uniprot by species name. Github does not accept file larger than 25MB so split into two files.
   5) `kegglist.txt`: KEGG prokaryote organism names and ID, for extracting the EC numbers from KEGG by species name.

scripts to run matlab locally
```bash
# extract unique metabolite names
while read line; do 
	echo "$line"
	arch -x86_64 /Applications/MATLAB_R2024b.app/bin/matlab \
      -nodisplay -nosplash -nodesktop < /Users/xinmengliao/Documents/Project/20241117_CMA_PD_ecGEMs/ecGEMs/${line}/meta-gecko-extractMet.m;
done < s.txt

# convert GEM to ecGEM
while read line; do 
	echo "$line"
	arch -x86_64 /Applications/MATLAB_R2024b.app/bin/matlab \
      -nodisplay -nosplash -nodesktop < /Users/xinmengliao/Documents/Project/20241117_CMA_PD_ecGEMs/ecGEMs/${line}/meta-gecko.m;
done < s.txt
```

### R-scripts 
1. Extracting species KEGG ID and Uniprot ID: `taxonomy.R`
2. Extracting SMILEs (PubChem API): `correctMetName.R`
3. Metagenomics analysis
4. Different abundance anlaysis for each species: `Differential abundance analysis.Rmd`
5. Constructing Metabolite-Microbe network and Reporter metabolite analysis: `ReporterMetabolites.Rmd`
