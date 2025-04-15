configfile: "GEMconfig.yaml"
HOME="/mnt/storage_pool/xinmeng/ecgem"

import os
import glob

def get_sampleids_from_path_pattern(path_pattern):
    ids = [os.path.basename(val).split('_R')[0] for val in (glob.glob(path_pattern))]
    new_list = []
    for var in ids:
        if var not in new_list:
            new_list.append(var)
    return new_list

SAMPLE_INDEX = get_sampleids_from_path_pattern('Study1/*')

rule all:
	input:
        	expand("{home}/prokka/{sampleID}", home=HOME, sampleID=SAMPLE_INDEX) ## Prokka


rule prokka:
	input:
		bins = '{home}/dna_bins/{sampleID}'
	output:
		directory('{home}/prokka/{sampleID}')
	shell:
		"""
        #conda activate prokka 
        #set +u;source activate {config[envs][prokkaroary]};set -u
        idvar=$(basename {output})
        mkdir -p {output}/
	mkdir -p {config[path][scratch]}/{config[folder][prokka]}/${{idvar}}/
	cp {input}/*.fa {config[path][scratch]}/{config[folder][prokka]}/${{idvar}}/
        cd {config[path][scratch]}/{config[folder][prokka]}/${{idvar}}/
        for bin in *.fa; do 
            echo ${{bin}}
            binID=$(echo ${{bin}} |sed "s/.fa//g")
            prokka -locustag ${{binID}} --cpus {config[cores][prokka]} --centre MAG --compliant -outdir ${{binID}} -prefix ${{binID}} ${{bin}} --force --quiet
            mv ${{binID}}/ {output};
            echo -e 'Prokka for ${{idvar}} ${{binID}} finished.'
        done

		"""
