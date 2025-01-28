#!/bin/bash
#These scripts were written by M Amine Hassani - hassani.medamine@gmail.com
#to analyze recombiation events in 82 B. burgdorferi sensu stricto with B31 ref genome
#module load anaconda
#source activate anvio-8
#source /home/mahassani/miniconda3/bin/activate anvio-8
### debug
set -euo pipefail
thread=4;
work_dir="/Users/mh057/Documents/BbSensuStricto/fastas/Gubbins"; # set working dir.
path_data="/Users/mh057/Documents/BbSensuStricto/fastas/Gubbins/fasta_files"; # path to fasta files
path_ref="/Users/mh057/Documents/BbSensuStricto/fastas/Gubbins/B31_ref/GCF_000008685.2_ASM868v2_genomic.fna";
#outfile=logfile
outfile=$work_dir"/log.txt"
errfile=$work_dir"/err.txt"
exec > >(cat >> $outfile)
exec 1> >(tee -a $outfile >&1)
exec 2> >(tee -a $errfile >&2)
log()
{
    echo "$1"
}

	cd $path_data

	log ">>> 1 reformat to *.fasta ..."
	cd $path_data;
	
	for i in fasta_files/*.fna; 
	do 
		name=$(basename "$i" .fasta);  echo -e "$name\t$PWD/$i"; 
	done > input.tab
	
	python $work_dir/generate_ska_alignment.py --reference $path_ref --input $path_data/input.tab --out $path_data/multi_ska_out.aln --threads 3 --k 29; 
	
	run_gubbins.py -p LD83 -c 5 -m 5 -v  -s SGCs_LD83_rooted.tre --sh-test --p-value 0.01\
	--extensive-search --filter-percentage 60 --fcondairst-tree-builder fasttree\
	--bootstrap 100  --tree-builder raxmlng --first-model JC --model GTRGAMMA;
	
	python extract_gubbins_clade_statistics.py --clades clades.txt --gff Bbss_out_gbbins.recombination_predictions.gff --snps Bbss_out_gbbins.branch_base_reconstruction.embl --tree Bbss_out_gbbins.final_tree.tre --out gubbins/out --print-trees --print-rec-lengths

echo " *** analysis done *** "