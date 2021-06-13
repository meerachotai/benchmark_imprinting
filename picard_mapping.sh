#!/usr/bin/env bash
# 
# strainA="cviA"
# strainB="cviB"
# outdir="june12_6pm"
# scripts_dir="/u/scratch/m/mchotai/rnaseq_simul/scripts_import"
# genome="$( pwd )/$outdir/cviA_genome.fa"
# annot="$( pwd )/$outdir/cviA_annot_picard.gff3"
# picard="/u/scratch/m/mchotai/imprinting/imprinting_analysis-master"
# ${scripts_dir}/picard_mapping.sh -A cviA -B cviB -g $genome -a $annot -x 2 -y 3 -M 95 -P 25 -D $picard -o $outdir -O ${outdir}/picard -d $scripts_dir

outdir=""
picard=""
scripts_dir=""

strainA=""
strainB=""

snps=""
genome=""
annot=""

mat_cutoff=95
pat_cutoff=25

rep=""
AxB_rep=""
BxA_rep=""

majority=""

paired=false
outprefix="out"

while getopts "A:B:a:g:r:o:M:P:m:D:d:s:f:x:y:O:p" opt; do
	case $opt in
		A)	strainA="$OPTARG"
			;;
		B)	strainB="$OPTARG"
			;;
		a)	annot="$OPTARG"
			;;
		g)	genome="$OPTARG"
			;;
		r)	rep="$OPTARG"
			;;
		o)	outdir="$OPTARG"
			;;
		M)	mat_cutoff="$OPTARG"
			;;
		P)	pat_cutoff="$OPTARG"
			;;
		m)	majority="$OPTARG"
			;;
		D)	picard="$OPTARG"
			;;
		d)	scripts_dir="$OPTARG"
			;;
		s)	snps="$OPTARG"
			;;
		f)	fastq_dir="$OPTARG"
			;;
		x)	AxB_rep="$OPTARG"
			;;
		y)	BxA_rep="$OPTARG"
			;;
		O)	outprefix="$OPTARG"
			;;
		p)	paired="true"
			;;
	esac
done

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

if [ ${#snps} == 0 ]; then
	rm $outdir/${strainA}_${strainB}_SNP.txt
	# a BED file containing SNPs between AxB and BxA, where the SNP is given as [A allele]>[B allele]
	for f in $outdir/per_chrom/*.refseq2simseq.SNP.vcf 
	do 
		sed -e '1,12d' $f > tmp.txt # trim top 12 lines
		sed -e '$s/$//' tmp.txt >> $outdir/${strainA}_${strainB}_SNP.txt
		rm tmp.txt
	done  
	cat $outdir/${strainA}_${strainB}_SNP.txt | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $4 ">" $5}' > $outdir/${strainA}_${strainB}_SNP_fin.txt
	snps="$outdir/${strainA}_${strainB}_SNP_fin.txt"
fi

if [ ${#fastq_dir} == 0 ]; then
	fastq_dir="$( pwd )/$outdir/reads_simul/${strainA}_${strainB}_"
fi

workdir=$( pwd )
outdir=${workdir}/${outdir}

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
${picard}/make_metagenome.py $snps $genome $map/${strainA}_${strainB}_meta --GTF $annot
mkdir $map/${strainA}_${strainB}_meta_STAR
STAR --runMode genomeGenerate --outFileNamePrefix "$map/${strainA}_${strainB}_meta_STAR/log" --genomeDir "$map/${strainA}_${strainB}_meta_STAR" --genomeFastaFiles $map/${strainA}_${strainB}_meta.fa --sjdbGTFfile $map/${strainA}_${strainB}_meta_metagtf.gtf --sjdbOverhang 49

count=1

printf "Calling imprinting...\n"
if [ "$paired" = "true" ]; then
	printf "Assuming that reciprocal crosses are paired by replicate\n"
	
	# mapping
	for i in $(seq 1 1 $rep)
	do
		cross="AxB_${i}"
		${picard}/rna_seq_map.sh -1 ${fastq_dir}${cross}.fq -g $map/${strainA}_${strainB}_meta_STAR -C $map/${strainA}_${strainB}_meta_metachrom.txt -o $map/${cross} -A $strainA -B $strainB -n ${cross} -a GATCGGAAGAGCGGTTCAG -3 -r
		cross="BxA_${i}"
		${picard}/rna_seq_map.sh -1 ${fastq_dir}${cross}.fq -g $map/${strainA}_${strainB}_meta_STAR -C $map/${strainA}_${strainB}_meta_metachrom.txt -o $map/${cross} -A $strainA -B $strainB -n ${cross} -a GATCGGAAGAGCGGTTCAG -3 -r
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
		
		${picard}/call_imprinting.sh -o $map/rep_${i}_imprinting -1 $AxB_bam -2  $BxA_bam -S $snps -G $annot -A $strainA -B $strainB -n rep_${i} -R 2 -I 2 -C 10 -M $mat_cutoff -P $pat_cutoff -c 10 -r #add -r if retrial
		
		cat $map/rep_${i}_${j}_imprinting/imprinting/rep_${i}_${j}_imprinting_filtered_MEGs.txt | awk -v var="$count" '{print $0 "\t"var }' >> ${outprefix}_all_MEGs.txt
		cat $map/rep_${i}_${j}_imprinting/imprinting/rep_${i}_${j}_imprinting_filtered_PEGs.txt | awk -v var="$count" '{print $0 "\t"var }' >> ${outprefix}_all_PEGs.txt
	done
	
	count=$rep
	
else
	printf "Assuming that reciprocal crosses are not paired by replicate, finding combinations...\n"
	
	# mapping
	for i in $(seq 1 1 $AxB_rep)
	do
		cross="AxB_${i}"
		${picard}/rna_seq_map.sh -1 ${fastq_dir}${cross}.fq -g $map/${strainA}_${strainB}_meta_STAR -C $map/${strainA}_${strainB}_meta_metachrom.txt -o $map/${cross} -A $strainA -B $strainB -n ${cross} -a GATCGGAAGAGCGGTTCAG -3 -r
	done
	
	for i in $(seq 1 1 $BxA_rep)
	do
		cross="BxA_${i}"
		${picard}/rna_seq_map.sh -1 ${fastq_dir}${cross}.fq -g $map/${strainA}_${strainB}_meta_STAR -C $map/${strainA}_${strainB}_meta_metachrom.txt -o $map/${cross} -A $strainA -B $strainB -n ${cross} -a GATCGGAAGAGCGGTTCAG -3 -r
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
			
			${picard}/call_imprinting.sh -o $map/rep_${i}_${j}_imprinting -1 $AxB_bam -2  $BxA_bam -S $snps -G $annot -A $strainA -B $strainB -n rep_${i}_${j} -R 2 -I 2 -C 10 -M $mat_cutoff -P $pat_cutoff -c 10 -r #add -r if retrial
			
			cat $map/rep_${i}_${j}_imprinting/imprinting/rep_${i}_${j}_imprinting_filtered_MEGs.txt | awk -v var="$count" '{print $0 "\t"var }' >> ${outprefix}_all_MEGs.txt
			cat $map/rep_${i}_${j}_imprinting/imprinting/rep_${i}_${j}_imprinting_filtered_PEGs.txt | awk -v var="$count" '{print $0 "\t"var }' >> ${outprefix}_all_PEGs.txt
			
			((count++))
		done
	done
fi

printf "Total reciprocal pairs inspected: ${count}\n"

if [ ${#majority} == 0 ]; then
	$scripts_dir/find_consensus.py -i ${outprefix}_all_MEGs.txt -t $count -o ${outprefix}_consensus_MEGs.txt
	$scripts_dir/find_consensus.py -i ${outprefix}_all_PEGs.txt -t $count -o ${outprefix}_consensus_PEGs.txt
else
	printf "Using majority voting of >= ${majority} for consensus calls\n"
	$scripts_dir/find_consensus.py -i ${outprefix}_all_MEGs.txt -m $majority -o ${outprefix}_consensus_MEGs.txt
	$scripts_dir/find_consensus.py -i ${outprefix}_all_PEGs.txt -m $majority -o ${outprefix}_consensus_PEGs.txt
fi

te=$(date +%s); echo "Done. Time elapsed: $( displaytime $(($te - $ts)) )"
