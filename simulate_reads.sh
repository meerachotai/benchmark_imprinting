# sample command:
# scripts_dir="/u/scratch/m/mchotai/rnaseq_simul/scripts_import"
# strainA="cviA"
# strainB="cviB"
# outdir="new_out"
# refA="$outdir/cviA_genome.fa"
# refB="$outdir/cviB_genome.fa" 
# 
# $scripts_dir/simulate_reads.sh -A $strainA -B $strainB -x $refA -y $refB -d $scripts_dir -s 5 -u 30 -m 10 -p 10 -r 50 -R 3 -M 95 -P 25 -o $outdir

# Required arguments ------------------------
scripts_dir=""
strainA=""
strainB=""
outdir=""
refA=""
refB=""

# Defaults ------------------------------
unbiased=30
meg=10
peg=10
seed=5
read_length=50
rep=3
megs_bias=95
pegs_bias=25

# made-up genome
# refA="${strainA}_Stranscripts.fa"
# annotA="${refA}_annot.gff3"
# 
# made-up genome + introduced mutations
# refB="${strainB}_genome.fa"
# annotB="${refB}_annot.gff3"

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


while getopts "A:B:x:y:d:m:p:u:s:r:R:M:P:o:" opt; do
	case $opt in
		A)	strainA="$OPTARG"
			;;
		B)	strainB="$OPTARG"
			;;
		x)	refA="$OPTARG"
			;;
		y)	refB="$OPTARG"
			;;
		d)	scripts_dir="$OPTARG"
			;;
		m)	meg="$OPTARG"
			;;
		p)	peg="$OPTARG"
			;;
		u)	unbiased="$OPTARG"
			;;
		s)	seed="$OPTARG"
			;;
		r)	read_length="$OPTARG"
			;;
		R)	rep="$OPTARG"
			;;
		M) 	megs_bias="$OPTARG"
			;;
		P)	pegs_bias="$OPTARG"
			;;
		o)	outdir="$OPTARG"
			;;
	esac
done

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
$scripts_dir/simulated_read-counts.R --seed $seed --disp med --n $unbiased --nMEG $meg --nPEG $peg --MEGbias $megs_bias --PEGbias $pegs_bias --rep $rep "${outdir}/count_simul"

# DO NOT CHANGE -----------------------------------
AxB="${outdir}/count_simul_AxB.txt"
BxA="${outdir}/count_simul_BxA.txt"

printf "Simulating reads...\n"
mkdir ${outdir}/reads_simul
${scripts_dir}/reads_simul.R -a $AxB -b $BxA -A $refA -B $refB -p I -r $read_length -s $seed -R $rep ${outdir}/reads_simul/simul

# remove spaces for mapping
for f in $(ls -v ${outdir}/reads_simul/simul_AxB_*) ; do echo "$(awk '{$1=$1};1' $f)" > $f ; done
for f in $(ls -v ${outdir}/reads_simul/simul_BxA_*) ; do echo "$(awk '{$1=$1};1' $f)" > $f ; done

te=$(date +%s); echo "Done. Time elapsed: $( displaytime $(($te - $ts)) )"
