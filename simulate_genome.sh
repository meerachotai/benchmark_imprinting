#!/usr/bin/env bash

# SIMULATING GENOMES with given similarity score

# dependencies: R, seqkit (https://bioinf.shenwei.me/seqkit/download/), 
# simuG (https://github.com/yjx1217/simuG)(note: you need to enter directory of simuG in -D)
# for now, this script does everything in the current directory so navigate to the right directory first (might change that later)
# need Biostrings library for R
# 
# sample command(s):
# 
# ref="/u/scratch/m/mchotai/rnaseq_simul/col_simul/Cvi.chr.all.v2.0.fasta"
# annot="/u/scratch/m/mchotai/rnaseq_simul/col_simul/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3"
# scripts_dir="/u/scratch/m/mchotai/rnaseq_simul/scripts_import"
# simuG="${scripts_dir}/simuG"
# 
# outdir="new_out"

# ${scripts_dir}/simulate_genome.sh -r $ref -a $annot -A cviA -B cviB -D $simuG -d $scripts_dir -S 90 -s 2 -m 1 -i 0 -e 2.5 -r 1 -n 50 -p 2 -t 10 -T 10 -W 2 -v 5 -o $outdir -x cviA_genome.fa -y cviB_genome.fa -X cviA_annot.gff3 -Y cviB_annot.gff3

# Required arguments ------------------------
scripts_dir=""
# directory=""
simuG=""

# original genome to start with, needs to be input by user
strainA=""
ref=""
annot=""

strainB=""
outdir=""

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

refA=""
annotA=""
refB=""
annotB=""

workdir=$( pwd )

# Flag options -----------------------------
skip1=false
transcript_error=false

while getopts "o:A:B:x:y:X:Y:r:a:D:d:S:s:i:e:m:r:n:v:p:t:T:W:2g" opt; do
	case $opt in
		o)	outdir="$OPTARG"
			;;
		A)	strainA="$OPTARG"
			;;
		B)	strainB="$OPTARG"
			;;
		x)	refA="$OPTARG"
			;;
		y)	refB="$OPTARG"
			;;
		X)	annotA="$OPTARG"
			;;
		Y)	annotB="$OPTARG"
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

outdir=${workdir}/${outdir}

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
printf "for reproducibility, using seed = ${seed} \n"
printf "working in directory: ${outdir} \n\n"

mkdir $outdir

#################
### functions ###
#################

make_annot() { 
	strain=$1
	ref=$2
	scripts_dir=$3
	out=$4
	outdir=$5
	
	samtools faidx ${ref}
	cut -f1-2 ${ref}.fai > ${outdir}/${strain}/${strain}_length+id.txt
	${scripts_dir}/make_annot.R ${outdir}/${strain}/${strain}_length+id.txt ${out}
}

transcript_maker() { 
	strain=$1
	ref=$2
	annot=$3
	outdir=$4
	
	# wget https://github.com/broadinstitute/picard/releases/download/2.24.2/picard.jar
	java -jar picard.jar NormalizeFasta --INPUT $ref --OUTPUT edited_$ref --LINE_LENGTH 50
	mv edited_$ref $ref
	gffread -w ${outdir}/${strain}/${strain}_transcripts.fa -g $ref $annot
}

# https://unix.stackexchange.com/questions/27013/displaying-seconds-as-days-hours-mins-seconds
displaytime () {
  local T=$1
  local D=$((T/60/60/24))
  local H=$((T/60/60%24))
  local M=$((T/60%60))
  local S=$((T%60))
  [[ $D > 0 ]] && printf '%d days ' $D
  [[ $H > 0 ]] && printf '%d hours ' $H
  [[ $M > 0 ]] && printf '%d minutes ' $M
  [[ $D > 0 || $H > 0 || $M > 0 ]] && printf 'and '
  printf '%d seconds\n' $S
}

#################
### pipeline ####
#################

# mkdir $outdir

time_start=$(date)	# time run was started
ts=$(date +%s)	# time run was started (in seconds)

echo "Run start on: $time_start"

# ------------------ step 1: make strainA ref and annot -----------------------------------------

# using original ref and annot to make a smaller genome
# produces a smaller genome fasta file called - $strainA_genome.fa (refA from now on)
if [ "$skip1" == "true" ]; then
	printf "\nSkipping step 1 (making strainA genome) by user request...\n"
else
	
	mkdir ${outdir}/$strainA
	
	# produces strainA_genome.fa - which is our new genome
	
	if [ "$transcript_error" == "false" ]; then
		gffread -w ${outdir}/${strainA}/${strainA}_transcripts.fa -g $ref $annot
	else
		transcript_maker $strainA $ref $annot $outdir
	fi
	
	# extract sequence ids
	grep '^>' ${outdir}/${strainA}/${strainA}_transcripts.fa | awk '{print $0}' | sed 's/^>//' > ${outdir}/${strainA}/${strainA}_seq_ids.txt

	# choose genes randomly
	${scripts_dir}/inv_transform_sampling.R ${outdir}/${strainA}/${strainA}_seq_ids.txt -n $total_n -o ${outdir}/${strainA}/${strainA}_genes.txt -s $seed

	# find sequences for sampled genes
	seqkit grep -n -f ${outdir}/${strainA}/${strainA}_genes.txt ${outdir}/${strainA}/${strainA}_transcripts.fa -o ${outdir}/$refA
	printf "\nsimulated FASTA file for strainA: ${refA}"
	# makes annotation for smaller genome
	
	make_annot $strainA ${outdir}/${refA} ${scripts_dir} ${outdir}/$annotA $outdir
	awk '{if(NR==1){print $0} else{print $0"strainA"}}' ${outdir}/$annotA > temp # add 'strainA' at the end for later mapping
	mv temp ${outdir}/$annotA 
	printf "\nsimulated annotation file for strainA: ${annotA}\n"
fi

# -------------------------- step 2: edit strainA to make strainB with a certain similarity score ------------------------------------

# add mutations (snps + indels) and make strain B reference sequence: strainB_genome.fa
# options: score snp_score indel_score (indel start) extend_score (indel extended / gap) match_score snp:indel-ratio total_n

$scripts_dir/edit_genome.sh -A $strainA -a ${outdir}/$refA -B $strainB -b ${outdir}/$refB -D $simuG -d $scripts_dir -S $score -s $snp_score -i $indel_score -e $extend_score -m $match_score -r $ratio -n $total_n -v $seed -p $pausetime -t $trials_simuG -W $window -T $trials_reject -o $outdir > ${outdir}/simul_${strainB}_log.txt
printf "\nsimulated FASTA file for strainB: ${refB}"

# use make_annot for strainB:
make_annot $strainB ${outdir}/${refB} $scripts_dir ${outdir}/${annotB} $outdir
awk '{if(NR==1){print $0} else{print $0"strainB"}}' ${outdir}/${annotB} > temp  # add 'strainB' at the end for later mapping
mv temp ${outdir}/${annotB}
printf "\nsimulated annotation file for strainB: ${annotB}\n"

te=$(date +%s); echo "Done. Time elapsed: $( displaytime $(($te - $ts)) )"
