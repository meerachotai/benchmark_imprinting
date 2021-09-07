#!/usr/bin/env bash

source $1
outprefix="${outdir}/picard"
mat_cutoff=${mat_cutoff_picard}
pat_cutoff=${pat_cutoff_picard}

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

###################
##### pipeline ####
###################

time_start=$(date)	# time run was started
ts=$(date +%s)	# time run was started (in seconds)
echo "Run start on: $time_start"

# make snp file for simulated SNPs, snps option must be left blank
# if [ ${#snps} == 0 ]; then
# 	# a BED file containing SNPs between AxB and BxA, where the SNP is given as [A allele]>[B allele]
# 	for f in $outdir/per_chrom/*.refseq2simseq.map.txt
# 	do 
# 		tail -n +2 $f | grep SNP | awk '{print $1 "\t" $7 "\t" $7+1 "\t" $5 ">" $10}' > $outdir/${strainA}_${strainB}_SNPs.txt
# 		snps="$outdir/${strainA}_${strainB}_SNPs.txt"
# 	done
# fi

# special options for simulated files, otherwise cannot leave these blank
if [ ${#snps} == 0 ]; then
	rm $outdir/${strainA}_${strainB}_SNP.txt
	tail -n +2 $outdir/${strainB}/${strainB}_ref2sim.txt | grep SNP | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $4 ">" $8}' > $outdir/${strainA}_${strainB}_SNP.txt
	snps="$outdir/${strainA}_${strainB}_SNP.txt"
fi

if [ ${#fastq_dir} == 0 ]; then
	fastq_dir="$( pwd )/$outdir/reads_simul/${strainA}_${strainB}_"
fi

workdir=$( pwd )
outdir=${workdir}/${outdir}/picard

mkdir -p $outdir

mkdir $outdir/picard_map
map="$outdir/picard_map"

printf "\nSummary of calls:\n"
printf "using reference genome file: ${genome}\n"
printf "using annotation files: ${annot} \n"
printf "outdirectory: ${outdir}\n"
printf "using SNP file: ${snps} \n"
printf "using FASTQ directory: ${fastq_dir}\n"
printf "using MEGs-cutoff: ${mat_cutoff}, PEGs-cutoff: ${pat_cutoff} for imprinting\n\n"

# remove header
# sed -e '1d' $annot > $map/${strainA}_annot.gff3
# annot="$map/${strainA}_annot.gff3"
# 
# # add transcript_id
# cat $annot| awk '{print $1 "\t" $2-1 "\t" $2 "\t" $4 ">" $5}' > $outdir/${strainA}_${strainB}_SNP.txt
# cat $annot | awk '{print $9}' | cut -d'=' -f 2 > ids.txt
# cat ids.txt | awk '{print "transcript_id " $1}' > transcript_ids.txt
# paste -d "; " $annot transcript_ids.txt > $map/${strainA}_annot_fin.gff3
# 
# sed -i 's/ID=/gene_id=/' $map/${strainA}_annot_fin.gff3
# 
# annot="$map/${strainA}_annot_fin.gff3"

printf "Making metagenome...\n"
${picard}/make_metagenome.py $snps $genome $map/${strainA}_${strainB} --GTF $annot
mkdir $map/${strainA}_${strainB}_meta_STAR
STAR --runMode genomeGenerate --outFileNamePrefix "$map/${strainA}_${strainB}_meta_STAR/log" --genomeDir "$map/${strainA}_${strainB}_meta_STAR" --genomeFastaFiles $map/${strainA}_${strainB}.fa --sjdbGTFfile $map/${strainA}_${strainB}_metagtf.gtf --sjdbOverhang 49

count=1

printf "Calling imprinting...\n"
if [ "$paired" = "true" ]; then
	printf "Assuming that reciprocal crosses are paired by replicate\n"
	
	if [ ${#rep} == 0 ]; then
		if [ ${#AxB_rep} ==  ${#BxA_rep} ]; then
			rep=$AxB_rep
		else
			printf "Number of AxB and BxA replicates given do not match, please re-check, or get combinations by setting config file option PAIRED_RECIPROCAL_CROSSES = FALSE\n"
		fi
	fi
	
	# mapping
	for i in $(seq 1 1 $rep)
	do
		if [ "$paired_end" == "true" ]; then
			# hisat2 -k 20 -S ${map}/${strainA}_${strainB}_${cross}_map.sam -x ${map}/concat_${strainA}_${strainB} -1 ${fastq_dir}${cross}_1.fq -2 ${fastq_dir}${cross}_2.fq
			cross="AxB_${i}"
			${picard}/rna_seq_map.sh -1 ${fastq_dir}${cross}_1.fq -2 ${fastq_dir}${cross}_2.fq -g $map/${strainA}_${strainB}_meta_STAR -C $map/${strainA}_${strainB}_metachrom.txt -o $map/${cross} -A $strainA -B $strainB -n ${cross} -a GATCGGAAGAGCGGTTCAG -3 -r # >> $map/mapping_log.txt
			cross="BxA_${i}"
			${picard}/rna_seq_map.sh -1 ${fastq_dir}${cross}_1.fq -2 ${fastq_dir}${cross}_2.fq -g $map/${strainA}_${strainB}_meta_STAR -C $map/${strainA}_${strainB}_metachrom.txt -o $map/${cross} -A $strainA -B $strainB -n ${cross} -a GATCGGAAGAGCGGTTCAG -3 -r # >> $map/mapping_log.txt
		else
			cross="AxB_${i}"
			${picard}/rna_seq_map.sh -1 ${fastq_dir}${cross}.fq -g $map/${strainA}_${strainB}_meta_STAR -C $map/${strainA}_${strainB}_metachrom.txt -o $map/${cross} -A $strainA -B $strainB -n ${cross} -a GATCGGAAGAGCGGTTCAG -3 -r # >> $map/mapping_log.txt
			cross="BxA_${i}"
			${picard}/rna_seq_map.sh -1 ${fastq_dir}${cross}.fq -g $map/${strainA}_${strainB}_meta_STAR -C $map/${strainA}_${strainB}_metachrom.txt -o $map/${cross} -A $strainA -B $strainB -n ${cross} -a GATCGGAAGAGCGGTTCAG -3 -r # >> $map/mapping_log.txt
		fi
	done
	
	# counting and calling
	rm ${outprefix}_all_MEGs.txt
	rm ${outprefix}_all_PEGs.txt
	for i in $(seq 1 1 $rep)
	do
		AxB="AxB_${i}"
		BxA="BxA_${i}"
		AxB_bam="${map}/${AxB}/STAR/${AxB}_unique_alignments.bam" 
		BxA_bam="${map}/${BxA}/STAR/${BxA}_unique_alignments.bam"
		
		${picard}/call_imprinting.sh -o $map/rep_${i}_${i}_imprinting -1 $AxB_bam -2  $BxA_bam -S $snps -G $annot -A $strainA -B $strainB -n rep_${i}_${i} -R 2 -I 2 -C 10 -M $mat_cutoff -P $pat_cutoff -c 10 -r >> $map/call_imprint_log.txt
		
		cat $map/rep_${i}_${i}_imprinting/imprinting/rep_${i}_${i}_imprinting_filtered_MEGs.txt | awk -v var="$count" '{print $0 "\t"var }' >> ${outprefix}_all_MEGs.txt
		cat $map/rep_${i}_${i}_imprinting/imprinting/rep_${i}_${i}_imprinting_filtered_PEGs.txt | awk -v var="$count" '{print $0 "\t"var }' >> ${outprefix}_all_PEGs.txt
	done
	
	count=$rep
	
else
	printf "Assuming that reciprocal crosses are not paired by replicate, finding combinations...\n"
	
	if [ ${#AxB_rep} == 0 ] || [ ${#BxA_rep} == 0 ]; then
		if [ ${#rep} == 0 ]; then
			printf "Number of replicates not given, please try again and use options -x and -y or -r"
		else
			AxB_rep=$rep
			BxA_rep=$rep
		fi
	fi
	
	echo Mapping...
	# mapping
	for i in $(seq 1 1 $AxB_rep)
	do
		cross="AxB_${i}"
		${picard}/rna_seq_map.sh -1 ${fastq_dir}${cross}.fq -g $map/${strainA}_${strainB}_meta_STAR -C $map/${strainA}_${strainB}_metachrom.txt -o $map/${cross} -A $strainA -B $strainB -n ${cross} -a GATCGGAAGAGCGGTTCAG -3 -r # >> $map/mapping_log.txt
	done
	
	for i in $(seq 1 1 $BxA_rep)
	do
		cross="BxA_${i}"
		${picard}/rna_seq_map.sh -1 ${fastq_dir}${cross}.fq -g $map/${strainA}_${strainB}_meta_STAR -C $map/${strainA}_${strainB}_metachrom.txt -o $map/${cross} -A $strainA -B $strainB -n ${cross} -a GATCGGAAGAGCGGTTCAG -3 -r # >> $map/mapping_log.txt
	done
	
	# counting and calling
	rm ${outprefix}_all_MEGs.txt
	rm ${outprefix}_all_PEGs.txt
	for i in $(seq 1 1 $AxB_rep); do
		for j in $(seq 1 1 $BxA_rep); do
			printf "Running combination ${count}: AxB - ${i} and BxA - ${j}\n"
			AxB="AxB_${i}"
			BxA="BxA_${j}"
			AxB_bam="${map}/${AxB}/STAR/${AxB}_unique_alignments.bam" 
			BxA_bam="${map}/${BxA}/STAR/${BxA}_unique_alignments.bam"
			
			${picard}/call_imprinting.sh -o $map/rep_${i}_${j}_imprinting -1 $AxB_bam -2  $BxA_bam -S $snps -G $annot -A $strainA -B $strainB -n rep_${i}_${j} -R 2 -I 2 -C 10 -M $mat_cutoff -P $pat_cutoff -c 10 -r # >> $map/call_imprint_log.txt
			
			cat $map/rep_${i}_${j}_imprinting/imprinting/rep_${i}_${j}_imprinting_filtered_MEGs.txt | awk -v var="$count" '{print $0 "\t"var }' >> ${outprefix}_all_MEGs.txt
			cat $map/rep_${i}_${j}_imprinting/imprinting/rep_${i}_${j}_imprinting_filtered_PEGs.txt | awk -v var="$count" '{print $0 "\t"var }' >> ${outprefix}_all_PEGs.txt
			
			((count++))
		done
	done
fi

printf "Total reciprocal pairs inspected: ${count}\n"

if [ ${#majority} == 0 ]; then
	$scripts_dir/find_consensus.py -i ${outprefix}_all_MEGs.txt -t $count -o ${outprefix}_MEGs.txt
	$scripts_dir/find_consensus.py -i ${outprefix}_all_PEGs.txt -t $count -o ${outprefix}_PEGs.txt
else
	printf "Using majority voting of >= ${majority} for consensus calls\n"
	$scripts_dir/find_consensus.py -i ${outprefix}_all_MEGs.txt -m $majority -o ${outprefix}_MEGs.txt
	$scripts_dir/find_consensus.py -i ${outprefix}_all_PEGs.txt -m $majority -o ${outprefix}_PEGs.txt
fi

te=$(date +%s); echo "Done. Time elapsed: $( displaytime $(($te - $ts)) )"
