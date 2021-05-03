#!/usr/bin/env bash

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
trials_reject=15

# Required arguments 
strainA=""
ref=""
strainB=""
simuG_dir=""
scripts_dir=""
total_n=""
outdir=$( pwd )


while getopts "A:a:B:b:D:d:S:s:i:e:m:r:n:v:p:t:T:W:o:l:" opt; do
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
		T)	trials_reject="$OPTARG"
			;;
		W)	window="$OPTARG"
			;;	
		o) 	outdir="$OPTARG"
			;;
		l)	low_range="$OPTARG"
			;;
	esac
done

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

time_start=$(date)	# time run was started
ts=$(date +%s)	# time run was started (in seconds)

echo "Run start on: $time_start"

per_chrom=$outdir/per_chrom
mkdir ${per_chrom}
mkdir ${outdir}/$strainB

rm ${per_chrom}/scores_log.txt

RANDOM=$seed # for reproducibility

printf "chromosomes unedited\n" > ${per_chrom}/un_edited_chr.txt
for i in $(seq $bottom_range 1 $total_n)
do
	grep '^>' ${refA} | awk '{print $0}' | sed 's/^>//' | sed "${i}q;d" > ${per_chrom}/chr_name.txt
	seqkit grep -n -f ${per_chrom}/chr_name.txt ${refA} > ${per_chrom}/${strainA}_${i}.fa

	$scripts_dir/scoring.R -A ${per_chrom}/${strainA}_${i}.fa -i $indel_score -e $extend_score -s $snp_score -m $match_score -p $score -r $ratio -P > ${per_chrom}/scoring_log_1.txt
	
	indel_count=$(grep indel-count: ${per_chrom}/scoring_log_1.txt | awk '{print $2}')
	snp_count=$(grep snp-count: ${per_chrom}/scoring_log_1.txt | awk '{print $2}')
	
	snp_log_score=$(grep snp-score-needed: ${per_chrom}/scoring_log_1.txt | awk '{print $2}')
	indel_log_score=$(grep indel-score-needed: ${per_chrom}/scoring_log_1.txt | awk '{print $2}')
	total_score=$(grep required-score: ${per_chrom}/scoring_log_1.txt | awk '{print $2}')
	
# 	printf "${indel_score}\n ${snp_score} \n ${total_score}\n"
	
	printf "chromosome ${i}:\n" >> ${per_chrom}/scores_log.txt
	printf "indel score needed: ${indel_log_score}\n" >> ${per_chrom}/scores_log.txt
	printf "snp score needed: ${snp_log_score}\n" >> ${per_chrom}/scores_log.txt
	printf "total score needed: ${total_score}\n" >> ${per_chrom}/scores_log.txt
	printf "trial\tindel score\tactual score\n" >> ${per_chrom}/scores_log.txt
	
	echo "chromosome ${i}" > ${per_chrom}/simuG_log.txt # refresh every time, so you can find the error for the particular chromosome
	counter_1=0
	stop=false
	echo try-again! > ${per_chrom}/scoring_log_2.txt # needed to begin the while loop below
	while grep try-again! ${per_chrom}/scoring_log_2.txt
	do
		if [ "$counter_1" = "$trials_reject" ]; then
			echo "Reached limit ${trials_reject} for rejection sampling, please revise limit / change window / change seed"
			printf "${i} \n" >> ${per_chrom}/un_edited_chr.txt # keeps track of chromosomes that reached a limit for rejection sampling
			break
		fi
		if [ "$stop" = "false" ]; then
			counter_2=1
			# Adapted from CLP's compare_BED_regions.sh code snippet:
			cmd="perl ${simuG_dir}/simuG.pl -r ${per_chrom}/${strainA}_${i}.fa -indel_count ${indel_count} -snp_count ${snp_count} -prefix ${per_chrom}/${strainB}_${i} -seed $RANDOM"
			timeout ${pausetime} $cmd >> ${per_chrom}/simuG_log.txt
			while [ $? -ne 0 ]; do
				if [ "$counter_2" = "$trials_simuG" ]; then
					stop=true
					break
				fi
				((counter_2++))
				echo "Retrying simuG... (trial ${counter_2++})"
				cmd="perl ${simuG_dir}/simuG.pl -r ${per_chrom}/${strainA}_${i}.fa -indel_count ${indel_count} -snp_count ${snp_count} -prefix ${per_chrom}/${strainB}_${i} -seed $RANDOM"
				timeout ${pausetime} $cmd >> ${per_chrom}/simuG_log.txt
			done
			cmd1="$scripts_dir/scoring.R -A ${per_chrom}/${strainA}_${i}.fa -i $indel_score -e $extend_score -s $snp_score -m $match_score -p $score -r $ratio -W $window -S ${per_chrom}/${strainB}_${i}.refseq2simseq.SNP.vcf -I ${per_chrom}/${strainB}_${i}.refseq2simseq.INDEL.vcf -V"
			$cmd1 > ${per_chrom}/scoring_log_2.txt
			
			((counter_1++))
			i_score=$(grep indel-score: ${per_chrom}/scoring_log_2.txt | awk '{print $2}')
			a_score=$(grep actual-score: ${per_chrom}/scoring_log_2.txt | awk '{print $2}')
# 			printf "${i_score} \n${a_score}\n"
			printf "${counter_1}\t${i_score}\t${a_score}\n" >> ${per_chrom}/scores_log.txt
			
		else
			echo "the similarity score = ${score}% is too low, please revise, refer per_chrom/scores_log.txt for details"
			break
		fi	
	done
done

for f in $(ls -v ${per_chrom}/*.simseq.genome.fa) ; do sed -e '$s/$//' $f ; done > ${refB}	

te=$(date +%s); echo "Done. Time elapsed: $( displaytime $(($te - $ts)) )"

# rm per_chrom/seed_reproducibility.txt
# echo "chosen seed: ${seed}" > per_chrom/seed_reproducibility.txt
# grep random\ seed: per_chrom/simuG_log.txt >> per_chrom/seed_reproducibility.txt	
