#!/usr/bin/env bash

source $1
outprefix="anderson"
mat_cutoff=${mat_cutoff_anderson}
pat_cutoff=${pat_cutoff_anderson}
pval=${pval_anderson}	
logfc=${logfc_anderson}
	 	
# ---------------------- step 7: rename annot and ref ------------------

# strain_name ref annot strain_type(A/B)
rename_chr() { 
	strain=$1
	ref=$2
	annot=$3
	map=$4
	
	# making new temp files, not changing original
	sed "s/^>/>${strain}_/g" $ref > $map/${strain}_genome.fa # any line that has > in it
	sed "s/^/${strain}_/" $annot > $map/${strain}_annot.gff3 # start of line	
	
	head="##gff-version 3" # replace header again
	sed -i "1s/.*/$head/" $map/${strain}_annot.gff3
}

# strainA strainB cross_name(AxB_rep/BxA_rep)
map() { 
	strainA=$1
	strainB=$2
	cross=$3
	map=$4	
	fastq_dir=$5
	paired_end=$6
	
	# concatenate reads files for A,B - moved to simulate_reads.sh
	# cat ${fastq_dir}${cross}_A.fq ${fastq_dir}${cross}_B.fq > ${map}/${strainA}_${strainB}_${cross}.fq
	if [ "$paired_end" == "true" ]; then
		hisat2 -k 20 -S ${map}/${strainA}_${strainB}_${cross}_map.sam -x ${map}/concat_${strainA}_${strainB} -1 ${fastq_dir}${cross}_1.fq -2 ${fastq_dir}${cross}_2.fq # --phred33
	else
		hisat2 -k 20 -S ${map}/${strainA}_${strainB}_${cross}_map.sam -x ${map}/concat_${strainA}_${strainB} ${fastq_dir}${cross}.fq # --phred33
	fi
}

# -------------------- step 10: counting -------------------------

# strainA strainB cross_name(AxB/BxA)
count() { 
	strainA=$1
	strainB=$2
	cross=$3
	map=$4
	htseq_i=$5
	
	python3 -m HTSeq.scripts.count -s no -m union -a 0 -i ${htseq_i} -o ${map}/${strainA}_${strainB}_${cross}_count.sam ${map}/${strainA}_${strainB}_${cross}_map.sam ${map}/concat_${strainA}_${strainB}.gff3  --nonunique all > ${map}/counts_${strainA}_${strainB}_${cross}.txt	
}

# # qsub -V -N counts -cwd -j y -o qsub_logs/counts.txt -m bae -b y -l h_rt=01:00:00,h_data=8G "$cmd"

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

time_start=$(date)	# time run was started
ts=$(date +%s)	# time run was started (in seconds)

echo "Run start on: $time_start"

printf "\nMake sure the following dependencies are loaded: hisat2, htseqcount\n"
printf "\nSummary of calls:\n"
printf "using reference genome files: ${refA}, ${refB}\n"
printf "using annotation files: ${annotA}, ${annotB}\n"
printf "outdirectory: ${outdir}\n"
printf "using FASTQ directory: ${fastq_dir}\n\n"
printf "Using HTSeq -i option as: ${htseq_i}"

workdir=$( pwd )
outdir=${workdir}/${outdir}

mkdir -p $outdir
mkdir ${outdir}/map
map="${outdir}/map" # where intermediate files will be stored

# preparing for concatenating:
# renames the ids to give strain names, adds A or B as needed for HTseq count

if [ "$rename" == "true" ]; then
	printf "Renaming chromosomes to match strainA/B...\n"
	rename_chr $strainA $refA $annotA $map
	refA="$map/${strainA}_genome.fa"
	annotA="$map/${strainA}_annot.gff3"

	rename_chr $strainB $refB $annotB $map
	refB="$map/${strainB}_genome.fa"
	annotB="$map/${strainB}_annot.gff3"
fi

# concatenates annot, ref + builds index for hisat2 (note: indexing takes long, possibly should qsub)
# NOTE: these annotations / refs need to be RENAMED with strainA and strainB names.
cat $refA $refB > $map/concat_${strainA}_${strainB}.fa
# { cat $refA; sed '1d' $refB; } > $map/concat_${strainA}_${strainB}.fa

hisat2-build $map/concat_${strainA}_${strainB}.fa $map/concat_${strainA}_${strainB} # build index
cat $annotA $annotB > $map/concat_${strainA}_${strainB}.gff3

# concatenates FASTQ files from reciprocal directions, maps them
# produces file strainA_strainB_cross.sam (for igv + HTseqcount)

# REQUIRED:
# give a 'fastq_dir' variable name that represents the start of the file name (including address to that file)
# the file HAS TO END with _$cross_$rep_A.fq or _$cross_$rep_B.fq, example: simul_AxB_1_A.fq

if [ ${#fastq_dir} == 0 ]; then
	fastq=()
	while IFS= read -r line; do
	  fastq+=("$line")
	done < $config
	printf "Mapping...\n"
	for i in $(seq 0 1 $(($rep - 1))); do
		if [ "$paired_end" == "true" ]; then
			start=$(($i * 4))
			cross=AxB_$(( $i + 1 ))
			hisat2 -k 20 -S ${map}/${strainA}_${strainB}_${cross}_map.sam -x ${map}/concat_${strainA}_${strainB} -1 ${fastq[${start}]} -2 ${fastq[$((${start} + 1))]} # --phred33
			cross=BxA_$(( $i + 1 ))
			hisat2 -k 20 -S ${map}/${strainA}_${strainB}_${cross}_map.sam -x ${map}/concat_${strainA}_${strainB} -1 ${fastq[$((${start} + 2))]} -2 ${fastq[$((${start} + 3))]} # --phred33
		else
			start=$(($i * 2))
			cross=AxB_$(( $i + 1 ))
			hisat2 -k 20 -S ${map}/${strainA}_${strainB}_${cross}_map.sam -x ${map}/concat_${strainA}_${strainB} ${fastq[${start}]} # --phred33
			cross=BxA_$(( $i + 1 ))
			hisat2 -k 20 -S ${map}/${strainA}_${strainB}_${cross}_map.sam -x ${map}/concat_${strainA}_${strainB} ${fastq[$((${start} + 1))]} # --phred33
		fi
	done
else
	printf "Mapping...\n"
	for i in $(seq 1 1 $rep)
	do
		map $strainA $strainB AxB_${i} $map $fastq_dir
		map $strainA $strainB BxA_${i} $map $fastq_dir
	done
fi

printf "Counting...\n"
# produces file strainA_strainB_cross_edit.sam (for igv)
# produces counts file strainA_strainB_cross.txt
for i in $(seq 1 1 $rep)
do
	count $strainA $strainB AxB_${i} $map $htseq_i
	count $strainA $strainB BxA_${i} $map $htseq_i
done

# gene key format MUST BE A | B, tab-delimited, syntelogs side-by-side, if user-provided
if [ ${#gene_key} == 0 ]; then
	printf "Creating gene key list...\n"
	printf "A\tB\n" > ${map}/${strainA}_${strainB}_gene_key.txt # for later
# 	paste <(awk 'BEGIN { FS="="} NR>1 { print $2 }' ${annotA}) <(awk 'BEGIN { FS="="} NR>1 { print $2 }' ${annotB}) >> ${map}/${strainA}_${strainB}_gene_key.txt
	paste <(cat $annotA | rev | cut -f 1 | rev | cut -d "=" -f 2-) <(cat $annotB | rev | cut -f 1 | rev | cut -d "=" -f 2-) >> ${map}/${strainA}_${strainB}_gene_key.txt
	sed -i -e "2d" ${map}/${strainA}_${strainB}_gene_key.txt
	gene_key="${map}/${strainA}_${strainB}_gene_key.txt"
fi

printf "Calling imprinting...\n"
# ${scripts_dir}/call_imprinting_anderson.R -c ${map}/counts_${strainA}_${strainB}_ -k $gene_key -p $pval -u $mat_cutoff -l $pat_cutoff -r $rep -f $logfc -A $strainA -B $strainB -a $a_annot -b $b_annot -C $outprefix

${scripts_dir}/get_counts_anderson.R -c ${map}/counts_${strainA}_${strainB}_ -k $gene_key -r $rep -A $strainA -B $strainB -a $a_annot -b $b_annot -C ${map}/${outprefix}
${scripts_dir}/call_imprinting_anderson.R -c ${map}/${outprefix}_ -p $pval -u $mat_cutoff -l $pat_cutoff -r $rep -f $logfc -A $strainA -B $strainB ${outdir}/${outprefix}

te=$(date +%s); echo "Done. Time elapsed: $( displaytime $(($te - $ts)) )"
