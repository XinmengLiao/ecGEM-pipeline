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
> For Linux 64, Open MPI is built with CUDA awareness but this support is disabled by default.                                                                 
To enable it, please set the environment variable OMPI_MCA_opal_cuda_support=true before                                                                     
launching your MPI processes. Equivalently, you can set the MCA parameter in the command line:                                                               
mpiexec --mca opal_cuda_support 1 ...   \                                                                                                                                                                                                                                                                               
> In addition, the UCX support is also built but disabled by default.                                                                                          
To enable it, first install UCX (conda install -c conda-forge ucx). Then, set the environment                                                                
variables OMPI_MCA_pml="ucx" OMPI_MCA_osc="ucx" before launching your MPI processes.                                                                         
Equivalently, you can set the MCA parameters in the command line:                                                                                            
mpiexec --mca pml ucx --mca osc ucx ...
> Note that you might also need to set UCX_MEMTYPE_CACHE=n for CUDA awareness via UCX.                                                                         
Please consult UCX's documentation for detail.
> 
> GTDB-Tk v1.7.0 requires ~40G of external data which needs to be downloaded                                                                              
     and unarchived. This can be done automatically, or manually:    
> 1. Run the command download-db.sh to automatically download to: /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem/share/gtdbtk-1.7.0/db/    \
> 2. Manually download the latest reference data:  https://github.com/Ecogenomics/GTDBTk#gtdb-tk-reference-data     \
> 3. (2b). Set the GTDBTK_DATA_PATH environment variable in the file: /cfs/klemming/projects/snic/naiss2023-23-637/ecgem/envs/metagem/etc/conda/activate.d


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
