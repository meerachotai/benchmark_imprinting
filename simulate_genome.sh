#!/usr/bin/env bash

# SIMULATING GENOMES with given similarity score

# dependencies: R, seqkit (https://bioinf.shenwei.me/seqkit/download/), 
# simuG (https://github.com/yjx1217/simuG)(note: you need to enter directory of simuG in -D)
# for now, this script does everything in the current directory so navigate to the right directory first (might change that later)
# 
# sample command(s):
# 
# ref="/u/scratch/m/mchotai/rnaseq_simul/col_simul/Cvi.chr.all.v2.0.fasta"
# annot="/u/scratch/m/mchotai/rnaseq_simul/col_simul/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3"
# scripts_dir="/u/scratch/m/mchotai/rnaseq_simul/scripts_import"
# simuG="${scripts_dir}/simuG"

# mkdir simul_trial
# cd simul_trial
# 
# ./simulate_genome_opt.sh -r $ref -a $annot -A cviA -B cviB -D $simuG -d $scripts_dir -S 70 -s 2 -m 1 -i 0 -e 2.5 -r 1 -n 50 -p 5 -t 10 -T 20 -W 1 -v 5

# Required arguments ------------------------
scripts_dir=""
# directory=""
simuG=""

# original genome to start with, needs to be input by user
strainA=""
ref=""
annot=""

strainB=""

# Defaults ------------------------------
seed=5

# based on BLAST summary
indel_score=0
extend_score=2.5
snp_score=2
match_score=1
ratio=1
seed=5
window=1
pausetime=5
trials_simuG=10
trials_reject=20
total_n=50

score=70 # to change, use -2

# Flag options -----------------------------
skip1=false
transcript_error=false

while getopts "A:B:r:a:D:d:S:s:i:e:m:r:n:v:p:t:T:W:2g" opt; do
	case $opt in
		A)	strainA="$OPTARG"
			;;
		B)	strainB="$OPTARG"
			;;
		r)	ref="$OPTARG"
			;;
		a)	annot="$OPTARG"
			;;
		D)	simuG="$OPTARG"
			;;
		d)	scripts_dir="$OPTARG"
			;;
		S)	score="$OPTARG"
			;;
		s)	snp_score="$OPTARG"
			;;
		i)	indel_score="$OPTARG"
			;;
		e)	extend_score="$OPTARG"
			;;
		m)	match_score="$OPTARG"
			;;
		r)	ratio="$OPTARG"
			;;
		n)	total_n="$OPTARG"
			;;
		v)	seed="$OPTARG"
			;;
		p)	pausetime="$OPTARG"
			;;
		t)	trials_simuG="$OPTARG"
			;;
		T)	trials_reject="$OPTARG"
			;;
		W)	window="$OPTARG"
			;;
		2)	# skip to step 2
			skip1=true
			;;
		g)	# error for gffread, redo with this flag
			transcript_error=true
			;;
	esac
done

printf "\nSummary of calls:\n"
printf "creating two genomes: ${strainA}, ${strainB}\n"
printf "Achieving a similarity score of ${score}%% \n"
printf "Using snp score = ${snp_score}, match score = ${match_score}, indel score: ${indel_score}, extend score: ${extend_score} \n"
printf "Simulating ${total_n} genes \n"
printf "simuG trial-kill settings: pausing for ${pausetime} seconds, for ${trials_simuG} times \n"
printf "for simulating strainB: rejection sampling for ${trials_reject} times and a window of ${window}%% \n"
if [ "$transcript_error" == "false" ]; then
	printf "running on mode: no error in getting transcripts from gffreads\n"
else
	printf "running on mode: error in getting transcripts from gffreads\n"
fi
printf "for reproducibility, using seed = ${seed} \n\n"


# DO NOT CHANGE -------------------------------
# made-up genome: strainA
refA="${strainA}_genome.fa"
annotA="${refA}_annot.gff3"

# made-up genome + introduced mutations: strainB
refB="${strainB}_genome.fa"
annotB="${strainB}_annot.gff3"

#################
### functions ###
#################

make_annot() { 
	ref=$1
	scripts_dir=$2

	samtools faidx ${ref}
	cut -f1-2 ${ref}.fai > ${ref}_length+id.txt
	${scripts_dir}/make_annot.R ${ref}_length+id.txt ${ref}
}

make_seq_table() { 
	strain=$1
	ref=$2

	# multi-line to single-line fasta
	awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ${strain}/${ref} > ${strain}/${strain}_sl_transcripts.fa
	# make table from .fa file
	cat ${strain}/${strain}_sl_transcripts.fa | paste - - > ${strain}/${strain}_seq.txt
}

transcript_maker() { 
	strain=$1
	ref=$2
	annot=$3
	
	# wget https://github.com/broadinstitute/picard/releases/download/2.24.2/picard.jar
	java -jar picard.jar NormalizeFasta --INPUT $ref --OUTPUT edited_$ref --LINE_LENGTH 50
	mv edited_$ref $ref
	gffread -w ${strain}/${strain}_transcripts.fa -g $ref $annot
}

#################
### pipeline ####
#################

# printf "\nMaking outdirectory ${directory}"

# mkdir ${directory}
# cd ${directory}

# ------------------ step 1: make strainA ref and annot -----------------------------------------

# using original ref and annot to make a smaller genome
# produces a smaller genome fasta file called - $strainA_genome.fa (refA from now on)
if [ "$skip1" == "true" ]; then
	printf "\nSkipping step 1 (making strainA genome) by user request...\n"
else
	
	mkdir $strainA
	
	# produces strainA_genome.fa - which is our new genome
	
	if [ "$transcript_error" == "false" ]; then
		gffread -w ${strainA}/${strainA}_transcripts.fa -g $ref $annot
	else
		transcript_maker $strainA $ref $annot
	fi
	
	# extract sequence ids
	grep '^>' ${strainA}/${strainA}_transcripts.fa | awk '{print $0}' | sed 's/^>//' > ${strainA}/${strainA}_seq_ids.txt

	# choose genes randomly
	${scripts_dir}/inv_transform_sampling.R ${strainA}/${strainA}_seq_ids.txt -n $total_n -o ${strainA}/${strainA}_genes.txt -s $seed

	# find sequences for sampled genes
	seqkit grep -n -f ${strainA}/${strainA}_genes.txt ${strainA}/${strainA}_transcripts.fa -o ${strainA}/${strainA}_genome.fa
	printf "\nsimulated FASTA file for strainA: ${refA}"
	# makes annotation for smaller genome
	
	make_annot ${strainA}/${refA} ${scripts_dir}
	printf "\nsimulated annotation file for strainA: ${annotA}\n"
	
	# uses ref files to make strainA_seq.txt tables (format: seq_id	fasta_seq)
	make_seq_table $strainA $refA
fi

# -------------------------- step 2: edit strainA to make strainB with a certain similarity score ------------------------------------

# add mutations (snps + indels) and make strain B reference sequence: strainB_genome.fa
# options: score snp_score indel_score (indel start) extend_score (indel extended / gap) match_score snp:indel-ratio total_n

$scripts_dir/edit_genome.sh -R $strainA -r $refA -N $strainB -D $simuG -d $scripts_dir -S $score -s $snp_score -i $indel_score -e $extend_score -m $match_score -r $ratio -n $total_n -v $seed -p $pausetime -t $trials_simuG -W $window -T $trials_reject > simul_${strainB}_log.txt
# qsub -V -N simulate_genome -cwd -j y -o qsub_logs/simul_genome.txt -m bae -b y -l h_rt=05:00:00,h_data=8G "$cmd"
printf "\nsimulated FASTA file for strainB: ${refB}"

# use make_annot for strainB:
make_annot $strainB/$refB $scripts_dir
printf "\nsimulated annotation file for strainB: ${annotB}\n"

# uses ref files to make strainB_seq.txt tables (format: seq_id	fasta_seq)
make_seq_table $strainB $refB

