default conda path and package path do not have enough space for installation. 

## Metagem conda env 
```
conda create -p /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem
```
```
conda activate /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem
```
```
vi /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem/.condarc
```
> pkgs_dirs:
>  - /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/pkgs
```
conda install metagem
```


## Metawrap conda env
```
conda create --prefix /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metawrap
```
```
conda activate /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metawrap
```
```
vi /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metawrap/.condarc
```
> pkgs_dirs:
>  - /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/pkgs
```
conda install -c ursky metawrap-mg=1.3.2
```
> Downloading and Extracting Packages:
                                                                                                                                                              
>Preparing transaction: done                                                                                                                                   
>Verifying transaction: done                                                                                                                                   
>Executing transaction: /                                                                                                                                      
>Krona installed.  You still need to manually update the taxonomy                                                                                              
>databases before Krona can generate taxonomic reports.  The update                                                                                            
>script is ktUpdateTaxonomy.sh.  The default location for storing                                                                                              
>taxonomic databases is /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metawrap/opt/krona/taxonomy                                                    
                                                                                                                                                              
>If you would like the taxonomic data stored elsewhere, simply replace                                                                                         
>this directory with a symlink.  For example:                                                                                                                  
                                                                                                                                                              
>rm -rf /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metawrap/opt/krona/taxonomy                                                                    
>mkdir /path/on/big/disk/taxonomy  
>ln -s /path/on/big/disk/taxonomy /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metawrap/opt/krona/taxonomy  
>ktUpdateTaxonomy.sh

>- The default QUAST package does not include:

>- GRIDSS (needed for structural variants detection)
>- SILVA 16S rRNA database (needed for reference genome detection in metagenomic datasets)
>- BUSCO tools and databases (needed for searching BUSCO genes) -- works in Linux only!

>To be able to use those, please run  
> quast-download-gridss  
> quast-download-silva  
> quast-download-busco

>done  
>ERROR: This cross-compiler package contains no program /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metawrap/bin/x86_64-conda_cos6-linux-gnu-gfortran  
>INFO: deactivate-gfortran_linux-64.sh made the following environmental changes:  
>-HOST=login1  
>ERROR: This cross-compiler package contains no program /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metawrap/bin/x86_64-conda_cos6-linux-gnu-gfortran  
>INFO: activate-gfortran_linux-64.sh made the following environmental changes:  
>+HOST=x86_64-conda_cos6-linux-gnu  
>-HOST=x86_64-conda-linux-gnu
