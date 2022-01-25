#!/usr/bin/env bash

source $1

logfc=${logfc_wyder}
pval=${pval_wyder}
outprefix="wyder"
counts_dir=${counts_dir_wyder}

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

printf "note: if you want to use your own counts, use call_imprinting_wyder.R directly\n"

echo "Run start on: $time_start"

printf "\nSummary of calls:\n"
printf "for ${strainA} and ${strainB}:\n"
printf "Wyder/edgeR settings: log2fc = ${logfc}, p-value = ${pval}\n\n"

workdir=$( pwd )
outdir=${workdir}/${outdir}

mkdir -p $outdir

mkdir $outdir/wyder
wyder="$outdir/wyder"

	
if [ ${#counts_dir} == 0 ]; then
	counts_dir=$outdir/picard/picard_map
fi

count=1
if [ "$paired" = "true" ]; then

	printf "Assuming that reciprocal crosses are paired by replicate, looking for count files in directory ${counts_dir}\n"
	
	if [ ${#rep} == 0 ]; then
		if [ ${#AxB_rep} ==  ${#BxA_rep} ]; then
			rep=$AxB_rep
		else
			printf "Number of AxB and BxA replicates given do not match, please re-check, or get combinations without using option -p\n"
		fi
	fi
	
	printf "Using replicates: ${rep}\n"
	for i in $(seq 1 1 $rep)
	do
		cp ${counts_dir}/rep_${i}_${i}_imprinting/counts_per_gene/rep_${i}_${i}_${strainA}x${strainB}_counts_merged.txt $wyder/wyder_AxB_${i}.txt
		cp ${counts_dir}/rep_${i}_${i}_imprinting/counts_per_gene/rep_${i}_${i}_${strainB}x${strainA}_counts_merged.txt $wyder/wyder_BxA_${i}.txt
	done
	count=$rep
	
else
	
	printf "Assuming that reciprocal crosses are not paired by replicate, looking for combinations count files in directory ${counts_dir}\n"
	
	if [ ${#AxB_rep} == 0 ] || [ ${#BxA_rep} == 0 ]; then
		if [ ${#rep} == 0 ]; then
			printf "Number of replicates not given, please try again and use options -x and -y or -r"
		else
			AxB_rep=$rep
			BxA_rep=$rep
		fi
	fi
	
	printf "Using AxB replicates: ${AxB_rep} and BxA replicates: ${BxA_rep}\n"
	# move over counts from Picard pipeline, rename to required format (_cross_rep.txt)
	for i in $(seq 1 1 $AxB_rep); do
		for j in $(seq 1 1 $BxA_rep); do
			cp $counts_dir/rep_${i}_${j}_imprinting/counts_per_gene/rep_${i}_${j}_${strainA}x${strainB}_counts_merged.txt $wyder/wyder_AxB_${count}.txt
			cp $counts_dir/rep_${i}_${j}_imprinting/counts_per_gene/rep_${i}_${j}_${strainB}x${strainA}_counts_merged.txt $wyder/wyder_BxA_${count}.txt
			((count++))
		done
	done
fi

((count--))
printf "Calling imprinted genes...\n"
$scripts_dir/call_imprinting_wyder.R -C -c "${wyder}/wyder_" -r $count -p $pval -l $logfc ${outdir}/$outprefix

te=$(date +%s); echo "Done. Time elapsed: $( displaytime $(($te - $ts)) )"

