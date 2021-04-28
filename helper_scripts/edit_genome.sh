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
ref_strain=""
ref=""
new_strain=""
simuG_dir=""
scripts_dir=""
total_n=""


while getopts "R:r:N:D:d:S:s:i:e:m:r:n:v:p:t:T:W:" opt; do
	case $opt in
		R)	ref_strain="$OPTARG"
			;;
		r)	ref="$OPTARG"
			;;
		N)	new_strain="$OPTARG"
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
	esac
done

mkdir per_chrom
mkdir $new_strain

rm per_chrom/indel_scores.txt

RANDOM=$seed # for reproducibility

printf "chromosomes unedited\n" > per_chrom/un_edited_chr.txt
for i in $(seq 1 1 $total_n)
do
	grep '^>' ${ref_strain}/${ref_strain}_genome.fa | awk '{print $0}' | sed 's/^>//' | sed "${i}q;d" > per_chrom/chr_name.txt
	seqkit grep -n -f per_chrom/chr_name.txt ${ref_strain}/${ref_strain}_genome.fa > per_chrom/${ref_strain}_${i}.fa

	$scripts_dir/scoring.R -A $ref_strain/${ref_strain}_seq.txt -i $indel_score -e $extend_score -s $snp_score -m $match_score -p $score -r $ratio -c $i -P > per_chrom/scoring_log_1.txt

	indel_count=$(grep indel-count: per_chrom/scoring_log_1.txt | awk '{print $2}')
	snp_count=$(grep snp-count: per_chrom/scoring_log_1.txt | awk '{print $2}')
	
	snp_score=$(grep snp-score-needed: per_chrom/scoring_log_1.txt | awk '{print $2}')
	indel_score=$(grep indel-score-needed: per_chrom/scoring_log_1.txt | awk '{print $2}')
	total_score=$(grep required-score: per_chrom/scoring_log_1.txt | awk '{print $2}')
	
	printf "chromosome ${i}: needed vs. actual comparison:\n" >> per_chrom/indel_scores.txt
	printf "indel score needed: ${indel_score}\n" >> per_chrom/indel_scores.txt
	printf "snp score needed: ${snp_score}\n" >> per_chrom/indel_scores.txt
	printf "total score needed: ${total_score}\n" >> per_chrom/indel_scores.txt
	printf "trial\tindel score\tactual score\n" >> per_chrom/indel_scores.txt
	
	echo "chromosome ${i}" > per_chrom/simuG_log.txt # refresh every time, so you can find the error for the particular chromosome
	counter_1=0
	stop=false
	echo try-again! > per_chrom/scoring_log_2.txt # needed to begin the while loop below
	while grep try-again! per_chrom/scoring_log_2.txt
	do
		if [ "$counter_1" = "$trials_reject" ]; then
			echo "Reached limit ${trials_reject} for rejection sampling, please revise limit / change window / change seed"
			printf "${i} \n" >> per_chrom/un_edited_chr.txt # keeps track of chromosomes that reached a limit for rejection sampling
			break
		fi
		if [ "$stop" = "false" ]; then
			counter_2=1
			# Adapted from CLP's compare_BED_regions.sh code snippet:
			cmd="perl ${simuG_dir}/simuG.pl -r per_chrom/${ref_strain}_${i}.fa -indel_count ${indel_count} -snp_count ${snp_count} -prefix per_chrom/${new_strain}_${i} -seed $RANDOM"
			timeout ${pausetime} $cmd >> per_chrom/simuG_log.txt
			while [ $? -ne 0 ]; do
				if [ "$counter_2" = "$trials_simuG" ]; then
					stop=true
					break
				fi
				((counter_2++))
				echo "Retrying simuG... (trial ${counter_2++})"
				cmd="perl ${simuG_dir}/simuG.pl -r per_chrom/${ref_strain}_${i}.fa -indel_count ${indel_count} -snp_count ${snp_count} -prefix per_chrom/${new_strain}_${i} -seed $RANDOM"
				timeout ${pausetime} $cmd >> per_chrom/simuG_log.txt
			done
			$scripts_dir/scoring.R -A $ref_strain/${ref_strain}_seq.txt -i $indel_score -e $extend_score -s $snp_score -m $match_score -p $score -r $ratio -c $i -W $window -S per_chrom/${new_strain}_${i}.refseq2simseq.SNP.vcf -I per_chrom/${new_strain}_${i}.refseq2simseq.INDEL.vcf -V > per_chrom/scoring_log_2.txt
			
			((counter_1++))
			i_score=$(grep indel-score: per_chrom/scoring_log_2.txt | awk '{print $2}')
			a_score=$(grep actual-score: per_chrom/scoring_log_2.txt | awk '{print $2}')
			printf "${counter_1}\t${i_score}\t${a_score}\n" >> per_chrom/indel_scores.txt
			
		else
			echo "the similarity score = ${score}% is too low, please revise"
			break
		fi	
	done
done

for f in per_chrom/*.simseq.genome.fa ; do sed -e '$s/$//' $f ; done > ${new_strain}/${new_strain}_genome.fa	

# rm per_chrom/seed_reproducibility.txt
# echo "chosen seed: ${seed}" > per_chrom/seed_reproducibility.txt
# grep random\ seed: per_chrom/simuG_log.txt >> per_chrom/seed_reproducibility.txt	
