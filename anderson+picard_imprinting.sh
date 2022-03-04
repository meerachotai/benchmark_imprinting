# param=counts
# source benchmark_files/benchmark_imprint_${param}.txt

source $1
pval=${pval_picard}
mat_cutoff=${mat_cutoff_picard}
pat_cutoff=${pat_cutoff_picard}

mkdir $outdir/map_picard
newMap=$outdir/map_picard
map=$outdir/map
count=$rep
outprefix=$outdir/anderson_picard

for i in $(seq 1 1 $rep); do
	AxB_A=$map/anderson_AxB_A_${i}.txt
	AxB_B=$map/anderson_AxB_B_${i}.txt
	BxA_A=$map/anderson_BxA_A_${i}.txt
	BxA_B=$map/anderson_BxA_B_${i}.txt

	${picard}/call_imprinting.sh -o $newMap/rep_${i}_${i}_imprinting -x $AxB_A -y $AxB_B -X $BxA_A -Y $BxA_B -A $strainA -B $strainB -n rep_${i}_${i} -R 2 -I 2 -C 10 -M $mat_cutoff -P $pat_cutoff -c 10 -r -p ${pval} >> $newMap/imprinting_log.txt
	
	cat $newMap/rep_${i}_${i}_imprinting/imprinting/rep_${i}_${i}_imprinting_filtered_MEGs.txt | awk -v var="$count" '{print $0 "\t"var }' >> ${outprefix}_all_MEGs.txt
	cat $newMap/rep_${i}_${i}_imprinting/imprinting/rep_${i}_${i}_imprinting_filtered_PEGs.txt | awk -v var="$count" '{print $0 "\t"var }' >> ${outprefix}_all_PEGs.txt
done



if [ ${#majority} == 0 ]; then
	$scripts_dir/find_consensus.py -i ${outprefix}_all_MEGs.txt -t $count -o ${outprefix}_MEGs.txt
	$scripts_dir/find_consensus.py -i ${outprefix}_all_PEGs.txt -t $count -o ${outprefix}_PEGs.txt
else
	printf "Using majority voting of >= ${majority} for consensus calls\n"
	$scripts_dir/find_consensus.py -i ${outprefix}_all_MEGs.txt -m $majority -o ${outprefix}_MEGs.txt
	$scripts_dir/find_consensus.py -i ${outprefix}_all_PEGs.txt -m $majority -o ${outprefix}_PEGs.txt
fi
