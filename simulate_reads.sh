# sample command:
# scripts_dir="/u/scratch/m/mchotai/rnaseq_simul/scripts_import"
# strainA="cviA"
# strainB="cviB"
# $scripts_dir/simulate_reads_opt.sh -A $strainA -B $strainB -d $scripts_dir -s 5 -u 35 -m 10 -p 15 -r 50 -R 3 -M 95 -P 25

# Required arguments ------------------------
scripts_dir=""
strainA=""
strainB=""

# Defaults ------------------------------
unbiased=30
meg=10
peg=10
seed=5
read_length=50
rep=3
megs_bias=95
pegs_bias=25

# DO NOT CHANGE -----------------------------------
AxB="count_simul_AxB.txt"
BxA="count_simul_BxA.txt"

# made-up genome
refA="${strainA}_Stranscripts.fa"
annotA="${refA}_annot.gff3"

# made-up genome + introduced mutations
refB="${strainB}_genome.fa"
annotB="${refB}_annot.gff3"

while getopts "A:B:d:m:p:u:s:r:R:M:P:" opt; do
	case $opt in
		A)	strainA="$OPTARG"
			;;
		B)	strainB="$OPTARG"
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
	esac
done

printf "\nSummary of calls:\n" 
printf "Simulating reads for genomes: ${strainA}, ${strainB}\n"
printf "Read length = ${read_length} \n"
printf "Number of replicates = ${rep} \n"
printf "for reproducibility, using seed = ${seed} \n\n"

printf "Simulating counts...\n"
$scripts_dir/simulated_read-counts.R --seed $seed --disp med --n $unbiased --nMEG $meg --nPEG $peg --MEGbias $megs_bias --PEGbias $pegs_bias --rep $rep "count_simul"

printf "Simulating reads...\n"
mkdir reads_simul
${scripts_dir}/reads_simul.R -a $AxB -b $BxA -A ${strainA}/${strainA}_seq.txt -B ${strainB}/${strainB}_seq.txt -p I -r $read_length -s $seed -R $rep reads_simul/simul


