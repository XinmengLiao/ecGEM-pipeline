configfile: "GEMconfig.yaml"
#SAMPLE_INDEX = {"Feces_20004_Visit3_S9"}
HOME="/cfs/klemming/projects/supr/naiss2023-23-637/ecgem/metaGEM/workflow"

import os
import glob

def get_sampleids_from_path_pattern(path_pattern):
    ids = [os.path.basename(val).split('_R')[0] for val in (glob.glob(path_pattern))]
    new_list = []
    for var in ids:
        if var not in new_list:
            new_list.append(var)
    return new_list

def get_binids_from_path_pattern(path_pattern):
    ids = [os.path.basename(val).split('.faa')[0] for val in (glob.glob(path_pattern))]
    new_list = []
    for var in ids:
        if var not in new_list:
            new_list.append(var)
    return new_list

#binIDs = get_binids_from_path_pattern('protein_bins/sample*/*.faa')
#GEMIDs = get_ids_from_path_pattern('GEMs/sample*/*.xml')
SAMPLE_INDEX = get_sampleids_from_path_pattern('dataset1/*')

rule all:
	input:
        	expand("{home}/qfiltered/{sampleID}/{sampleID}_R1_qfiltered.fastq.gz", home = HOME, sampleID = SAMPLE_INDEX), ## Fastp
        	expand("{home}/qfiltered/{sampleID}/{sampleID}_R2_qfiltered.fastq.gz", home = HOME, sampleID = SAMPLE_INDEX),
        	expand("{home}/megahit_assembled/{sampleID}/contigs.fasta.gz", home = HOME, sampleID = SAMPLE_INDEX) ## Megahit assemble



####--------------------------------------####
#### 1. Fastp (no problem, run smoothly)
####--------------------------------------####

rule fastp:
	input:
		R1 = "{home}/dataset1/{sampleID}/{sampleID}_R1.fastq.gz",
		R2 = "{home}/dataset1/{sampleID}/{sampleID}_R2.fastq.gz"
	output:
		R1 = "{home}/qfiltered/{sampleID}/{sampleID}_R1_qfiltered.fastq.gz",
		R2 = "{home}/qfiltered/{sampleID}/{sampleID}_R2_qfiltered.fastq.gz"
	threads: 32
	benchmark:
		"{home}/benchmark/{sampleID}.fastp.benchmark.txt"
	shell:
		"""
        echo -e "$(date)\nSection starts\n ***** Fastp ***** \n"

        #module load bioinfo-tools
        #module load conda
        #metagem # activate conda metagem env

        mkdir -p $(dirname {output.R1})
        fastp -i {input.R1} \
            -I {input.R2} \
            -o {output.R1} \
            -O {output.R2} \
	    --thread {config[cores][fastp]} \
            -j $(dirname {output.R1})/$(echo $(basename $(dirname {output.R1}))).json \
            -h $(dirname {output.R1})/$(echo $(basename $(dirname {output.R1}))).html

        echo "Fastp quality filtering done at $(date). "
		"""

 ## Skip -- qfilterVis


####--------------------------------------####
#### 2. Megahit (no problem, run smoothly)
####--------------------------------------####
# If end in the middle, can not directly restart
# how can set conda in general snakemake rule 

rule megahit:
    input:
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2
    output:
        "{home}/megahit_assembled/{sampleID}/contigs.fasta.gz"
    threads: 32
    benchmark:
        "{home}/benchmark/{sampleID}.megahit.benchmark.txt"
    log:
        "{home}/logs/{sampleID}/{sampleID}_megahit.log"
    shell:
        """
        echo -e "$(date)\nSection starts\n Megahit \n"
        
        # Actiavte the conda environment 
        #set +e
        #set +u 
        #source activate metagem 
        #set -u 
        #set -e

        mkdir -p $(dirname {output})

        #Create the temporary folder
        idvar=$(echo $(basename $(dirname {output})))
        ihome=$(echo $(dirname $(dirname {output})))
        echo -e "\nCreating temporary directory ${{ihome}}/${{idvar}} ..."
        mkdir -p ${{ihome}}/${{idvar}}
        cd ${{ihome}}/${{idvar}}

        echo -n "Copying qfiltered reads to ${{ihome}}/${{idvar}} ... "
        cp {input.R1} {input.R2} ./
        echo -n "Running megahit ..., assemble mode is {config[params][assemblyPreset]} ... "

        # make sure there is not previous intermediate file 
        rm -rf tmp 

        megahit -1 $(basename {input.R1}) \
            -2 $(basename {input.R2}) \
            -o tmp \
            -t {config[cores][megahit]} \
            --presets {config[params][assemblyPreset]} \
            --verbose \
            --min-contig-len {config[params][assemblyMin]} \
            #--k-min  \
            #--k-max 99 \
            #--k-step 20 \
            2> {log}
        echo "done. "
        echo "Renaming assembly ... "
        mv tmp/final.contigs.fa contigs.fasta

        # Remove spaces from the contig headers and replace with hyphens
        sed -i 's/ /-/g' contigs.fasta
        gzip contigs.fasta
	
	# remove tmp/ after completion
	rm -rf tmp/
        
        echo "Megahit assembling Done at $(date). "
	"""
