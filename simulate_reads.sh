#!/usr/bin/env bash

source $1

peg=$npeg
meg=$nmeg

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

time_start=$(date)	# time run was started
ts=$(date +%s)	# time run was started (in seconds)

echo "Run start on: $time_start"

workdir=$( pwd )
outdir=${workdir}/${outdir}

printf "\nSummary of calls:\n" 
printf "Simulating reads for genomes: ${strainA}, ${strainB}\n"
printf "using reference genome files: ${refA}, ${refB}\n"
printf "Read length = ${read_length} \n"
printf "Number of replicates = ${rep} \n"
printf "for reproducibility, using seed = ${seed} \n"
printf "outdirectory: ${outdir}\n\n"

printf "Simulating counts...\n"
$scripts_dir/simulate_counts.R --seed $seed --disp $disp --n $unbiased --nMEG $meg --nPEG $peg --MEGbias $megs_bias --PEGbias $pegs_bias --rep $rep --seqDepth $seq_depth "${outdir}/count_simul"

# DO NOT CHANGE -----------------------------------
AxB="${outdir}/count_simul_AxB.txt"
BxA="${outdir}/count_simul_BxA.txt"

printf "Simulating reads...\n"
mkdir ${outdir}/reads_simul
${scripts_dir}/reads_simul.R -a $AxB -b $BxA -A $outdir/${refA}.fa -B $outdir/${refB}.fa -p I -r $read_length -s $seed -R $rep ${outdir}/reads_simul/simul

# remove spaces for mapping
for f in $(ls -v ${outdir}/reads_simul/simul_AxB_*) ; do echo "$(awk '{$1=$1};1' $f)" > $f ; done
for f in $(ls -v ${outdir}/reads_simul/simul_BxA_*) ; do echo "$(awk '{$1=$1};1' $f)" > $f ; done

# concatenating both A and B reads to one file
for i in $(seq 1 1 $rep)
do
	cat ${outdir}/reads_simul/simul_AxB_${i}_A.fq ${outdir}/reads_simul/simul_AxB_${i}_B.fq > ${outdir}/reads_simul/${strainA}_${strainB}_AxB_${i}.fq
	cat ${outdir}/reads_simul/simul_BxA_${i}_A.fq ${outdir}/reads_simul/simul_BxA_${i}_B.fq > ${outdir}/reads_simul/${strainA}_${strainB}_BxA_${i}.fq
done

echo Preparing true MEGs and PEGs files...
while IFS= read -r line; do
sed -n "${line}p" ${outdir}/all_genes.txt >> ${outdir}/true_PEGs.txt 
done < ${outdir}/count_simul_pegs.txt
# sed -i 's/^.//' ${outdir}/true_pegs.txt
rm ${outdir}/count_simul_pegs.txt

while IFS= read -r line; do
sed -n "${line}p" ${outdir}/all_genes.txt >> ${outdir}/true_MEGs.txt 
done < ${outdir}/count_simul_megs.txt
# sed -i 's/^.//' ${outdir}/true_megs.txt
rm ${outdir}/count_simul_megs.txt

te=$(date +%s); echo "Done. Time elapsed: $( displaytime $(($te - $ts)) )"
