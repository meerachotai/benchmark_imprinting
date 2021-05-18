#!/usr/bin/env bash

# scripts_dir="/u/scratch/m/mchotai/rnaseq_simul/scripts_import"
# simuG_dir="${scripts_dir}/simuG"
# $scripts_dir/edit_genome_v2.sh -A cviA -a $( pwd )/new_out/cviA_genome.fa -B cviB -b alt_edit_genome/cviB_genome_alt.fa -D $simuG_dir -d $scripts_dir -S 90 -n 50 -v 5 -p 3 -t 10 -W 3 -o new_out

# Defaults
indel_score=0
extend_score=2.5
snp_score=2
match_score=1
ratio=1

score=70
window=1

seed=5
pausetime=5
trials_simuG=5

# Required arguments 
strainA=""
ref=""
strainB=""
simuG_dir=""
scripts_dir=""
total_n=""
outdir=$( pwd )

while getopts "A:a:B:b:D:d:S:s:i:e:m:r:n:v:p:t:W:o:l:" opt; do
	case $opt in
		A)	strainA="$OPTARG"
			;;
		a)	refA="$OPTARG"
			;;
		B)	strainB="$OPTARG"
			;;
		b)	refB="$OPTARG"
			;;
		D)	simuG_dir="$OPTARG"
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
		W)	window="$OPTARG"
			;;	
		o) 	outdir="$OPTARG"
			;;
		l)	low_range="$OPTARG"
			;;
	esac
done

RANDOM=$seed

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

printf "\nSummary of calls:\n"
printf "editing genome ${refA} to make genome for ${strainB}\n"
printf "Achieving a similarity score of ${score}%% \n"
printf "Using snp score = ${snp_score}, match score = ${match_score}, indel score: ${indel_score}, extend score: ${extend_score} \n"
printf "Simulating ${total_n} genes \n"
printf "simuG trial-kill settings: pausing for ${pausetime} seconds, for ${trials_simuG} times \n"
printf "for simulating strainB: rejection sampling with a window of ${window}%% \n"
printf "for reproducibility, using seed = ${seed} \n"
printf "working in directory: ${outdir} \n\n"

time_start=$(date)	# time run was started
ts=$(date +%s)	# time run was started (in seconds)

echo "Run start on: $time_start"

per_chrom=$outdir/per_chrom
mkdir ${per_chrom}
mkdir ${outdir}/$strainB

$scripts_dir/simuG_score.R -A ${refA} -i $indel_score -e $extend_score -s $snp_score -m $match_score -p $score -r $ratio -o ${per_chrom} > ${per_chrom}/simuG_score.txt

indel_count=$(grep indel-count: ${per_chrom}/simuG_score.txt | awk '{print $2}')
# printf $indel_count
snp_count=$(grep snp-count: ${per_chrom}/simuG_score.txt | awk '{print $2}')
# printf $snp_count

# since the R script appends to these files, 'refresh' them just in case
> ${per_chrom}/chromosomes_done.txt
> ${per_chrom}/score_stats.txt
# touch ${per_chrom}/chromosomes_done.txt
# req_len = $(wc -l < ${per_chrom}/per_chrom_stats.txt)

while [ $(wc -l < ${per_chrom}/chromosomes_done.txt) != $total_n ]; do
	cmd="perl ${simuG_dir}/simuG.pl -r ${refA} -indel_count ${indel_count} -snp_count ${snp_count} -prefix ${per_chrom}/${strainB} -seed $RANDOM" # sometimes defaults to 0, for some reason
# 	echo $cmd
	timeout ${pausetime} $cmd >> ${per_chrom}/simuG_log.txt
	counter=1
	while [ $? -ne 0 ]; do
		if [ "$counter" = "$trials_simuG" ]; then
			stop=true
			break
		fi
		((counter++))
		echo "Retrying simuG... (trial ${counter})"
		cmd="perl ${simuG_dir}/simuG.pl -r ${refA} -indel_count ${indel_count} -snp_count ${snp_count} -prefix ${per_chrom}/${strainB} -seed $RANDOM"
		timeout ${pausetime} $cmd >> ${per_chrom}/simuG_log.txt
	done
	$scripts_dir/vcf_score.R -G ${per_chrom}/${strainB}.simseq.genome.fa -B ${refB} -i $indel_score -e $extend_score -s $snp_score -W $window -S ${per_chrom}/${strainB}.refseq2simseq.SNP.vcf -I ${per_chrom}/${strainB}.refseq2simseq.INDEL.vcf -o ${per_chrom}
done

te=$(date +%s); echo "Done. Time elapsed: $( displaytime $(($te - $ts)) )"
