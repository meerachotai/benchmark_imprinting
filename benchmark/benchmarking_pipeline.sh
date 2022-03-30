#!/usr/bin/env bash

# Defaults: -------------

disp="med"

# FOR IMPRINTING
mat_bias=85
pat_bias=25
# FOR SIMULATING
meg_bias=90
peg_bias=20

nmeg=200
npeg=200
n=1000
seed=5
alpha=0.05
similarity=90
logfc=0
seqDepth=1

# Flag options --------------

falpha=false
fdisp=false
fmat=false
fpat=false
fsim=false
fseq=false

while getopts "o:mpadsD" opt; do
	case $opt in
		o)	outdir="$OPTARG"
			;;
		m)	# maternal-bias
			fmat=true
			;;
		p)	# paternal-bias
			fpat=true
			;;
		a)	# alpha
			falpha=true
			;;
		d)	# dispersion
			fdisp=true
			;;
		s)	# similarity
			fsim=true
			;;
		D)	# sequencing depth
			fseq=true
			;;
	esac
done

mkdir benchmark_files
mkdir benchmark_files/log

echo "outdir: ${outdir}"

################
### function ###
################

# runs ONLY IMPRINTING step
run_imprint() {
	parameter=$1
	j=$2
	
	echo for $parameter, using $j	

	benchmark_imprinting/anderson_imprinting.sh benchmark_files/benchmark_imprint_${parameter}.txt > benchmark_files/log/anderson_log_${parameter}.txt
	benchmark_imprinting/picard_imprinting.sh benchmark_files/benchmark_imprint_${parameter}.txt > benchmark_files/log/picard_log_${parameter}.txt
	benchmark_imprinting/wyder_imprinting.sh benchmark_files/benchmark_imprint_${parameter}.txt > benchmark_files/log/wyder_log_${parameter}.txt
	benchmark_imprinting/anderson+picard_imprinting.sh benchmark_files/benchmark_imprint_${param}.txt > benchmark_files/log/anderson+picard_log_${parameter}.txt
	benchmark_imprinting/roth_imprinting.sh benchmark_files/benchmark_imprint_${parameter}.txt > benchmark_files/log/roth_log_${parameter}.txt
	benchmark_imprinting/roth+wyder_imprinting.sh benchmark_files/benchmark_imprint_${parameter}.txt > benchmark_files/log/roth+wyder_log_${parameter}.txt
	benchmark_imprinting/anderson+wyder_imprinting.sh benchmark_files/benchmark_imprint_${parameter}.txt > benchmark_files/log/anderson+wyder_log_${parameter}.txt

	# for simulation, remove 'strainA'
	awk '{ gsub("strainA", "", $1); print }' $outdir/anderson_MEGs.txt > $outdir/tmp && mv $outdir/tmp $outdir/anderson_MEGs.txt
	awk '{ gsub("strainA", "", $1); print }' $outdir/anderson_PEGs.txt > $outdir/tmp && mv $outdir/tmp $outdir/anderson_PEGs.txt
	
	awk '{ gsub("strainA", "", $1); print }' $outdir/anderson_picard_MEGs.txt > $outdir/tmp && mv $outdir/tmp $outdir/anderson_picard_MEGs.txt
	awk '{ gsub("strainA", "", $1); print }' $outdir/anderson_picard_PEGs.txt > $outdir/tmp && mv $outdir/tmp $outdir/anderson_picard_PEGs.txt

	awk '{ gsub("strainA", "", $1); print }' $outdir/anderson_wyder_MEGs.txt > $outdir/tmp && mv $outdir/tmp $outdir/anderson_wyder_MEGs.txt
	awk '{ gsub("strainA", "", $1); print }' $outdir/anderson_wyder_PEGs.txt > $outdir/tmp && mv $outdir/tmp $outdir/anderson_wyder_PEGs.txt

	# output true positive stats
	benchmark_imprinting/benchmark/call_true_positives.R $outdir > benchmark_files/true_pos_${parameter}.txt
	
	meg=$(grep 'total\ Anderson\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_mat=$(grep 'true\ Anderson\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	peg=$(grep 'total\ Anderson\ paternally-biased:'  benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_pat=$(grep 'true\ Anderson\ paternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	
	printf "$j\t$meg\t$peg\t$true_mat\t$true_pat\n" >> benchmark_files/$parameter/anderson_${parameter}.txt

	meg=$(grep 'total\ Picard\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_mat=$(grep 'true\ Picard\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	peg=$(grep 'total\ Picard\ paternally-biased:'  benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_pat=$(grep 'true\ Picard\ paternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')

	printf "$j\t$meg\t$peg\t$true_mat\t$true_pat\n" >> benchmark_files/$parameter/picard_${parameter}.txt
	
	meg=$(grep 'total\ Anderson+Picard\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_mat=$(grep 'true\ Anderson+Picard\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	peg=$(grep 'total\ Anderson+Picard\ paternally-biased:'  benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_pat=$(grep 'true\ Anderson+Picard\ paternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')

	printf "$j\t$meg\t$peg\t$true_mat\t$true_pat\n" >> benchmark_files/$parameter/anderson_picard_combined_${parameter}.txt
	
	meg=$(grep 'total\ Anderson+Wyder\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_mat=$(grep 'true\ Anderson+Wyder\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	peg=$(grep 'total\ Anderson+Wyder\ paternally-biased:'  benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_pat=$(grep 'true\ Anderson+Wyder\ paternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')

	printf "$j\t$meg\t$peg\t$true_mat\t$true_pat\n" >> benchmark_files/$parameter/anderson_wyder_combined_${parameter}.txt

	meg=$(grep 'total\ Wyder\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_mat=$(grep 'true\ Wyder\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	peg=$(grep 'total\ Wyder\ paternally-biased:'  benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_pat=$(grep 'true\ Wyder\ paternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')

	printf "$j\t$meg\t$peg\t$true_mat\t$true_pat\n" >> benchmark_files/$parameter/wyder_${parameter}.txt
	
	meg=$(grep 'total\ Roth\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_mat=$(grep 'true\ Roth\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	peg=$(grep 'total\ Roth\ paternally-biased:'  benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_pat=$(grep 'true\ Roth\ paternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	
	printf "$j\t$meg\t$peg\t$true_mat\t$true_pat\n" >> benchmark_files/$parameter/roth_${parameter}.txt
	
	meg=$(grep 'total\ Roth+Wyder\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_mat=$(grep 'true\ Roth+Wyder\ maternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	peg=$(grep 'total\ Roth+Wyder\ paternally-biased:'  benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	true_pat=$(grep 'true\ Roth+Wyder\ paternally-biased:' benchmark_files/true_pos_${parameter}.txt | awk '{print $4}')
	
	printf "$j\t$meg\t$peg\t$true_mat\t$true_pat\n" >> benchmark_files/$parameter/roth_wyder_combined_${parameter}.txt
}

if [ "$fmat" == "true" ]; then

	echo "Changing %maternal for maternally-biased (for simulating)"
	param="matbias"
	array=(75 80 85 95 100 90)
# 	param="matbias_90"
# 	array=(90)
	mkdir benchmark_files/$param

	printf "Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/wyder_$param.txt
	printf "Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/picard_$param.txt
	printf "Anderson method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/anderson_$param.txt
	printf "Anderson+Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/anderson_picard_combined_$param.txt
	printf "Anderson+Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/anderson_wyder_combined_$param.txt
	printf "Roth method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/roth_$param.txt
	printf "Roth+Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/roth_wyder_combined_$param.txt

	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/picard_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/wyder_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_picard_combined_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_wyder_combined_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/roth_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/roth_wyder_combined_${param}.txt
	
	benchmark_imprinting/config/read_config_simul.py benchmark_imprinting/config/simulation_config.txt benchmark_files/benchmark_simul_${param}.txt

	benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_files/benchmark_imprint_${param}.txt
	benchmark_imprinting/benchmark/config_changer.sh -b $mat_bias -B $pat_bias -a $alpha -o $outdir -l $logfc -C benchmark_files/benchmark_imprint_${param}.txt
	
	for i in "${array[@]}"; do
	
		benchmark_imprinting/benchmark/config_changer.sh -n $n -m $i -p $peg_bias -M $nmeg -P $npeg -d $disp -s $similarity -o $outdir -D $seqDepth -c benchmark_files/benchmark_simul_${param}.txt

		benchmark_imprinting/simulate_genome.sh benchmark_files/benchmark_simul_${param}.txt overwrite #> benchmark_files/log/genome_log.txt
		benchmark_imprinting/simulate_reads.sh benchmark_files/benchmark_simul_${param}.txt #> benchmark_files/log/reads_log.txt
	
		# call imprinting
		run_imprint $param $i
	done
fi

if [ "$fpat" == "true" ]; then
	echo "Changing %maternal for paternally-biased (for simulating)"
	param="patbias"
	array=(0 10 20 30 40 45 50)
# 	param="patbias_0"
# 	array=(0)
	mkdir benchmark_files/$param
	
	printf "Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/wyder_$param.txt
	printf "Anderson method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/anderson_$param.txt
	printf "Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/picard_$param.txt
	printf "Anderson+Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/anderson_picard_combined_$param.txt
	printf "Anderson+Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/anderson_wyder_combined_$param.txt
	printf "Roth method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/roth_$param.txt
	printf "Roth+Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/roth_wyder_combined_$param.txt

	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/picard_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/wyder_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_picard_combined_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_wyder_combined_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/roth_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/roth_wyder_combined_${param}.txt

	benchmark_imprinting/config/read_config_simul.py benchmark_imprinting/config/simulation_config.txt benchmark_files/benchmark_simul_${param}.txt

	benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_files/benchmark_imprint_${param}.txt
	benchmark_imprinting/benchmark/config_changer.sh -b $mat_bias -B $pat_bias -a $alpha -o $outdir -l $logfc -C benchmark_files/benchmark_imprint_${param}.txt
	
	for i in "${array[@]}"; do
	
		benchmark_imprinting/benchmark/config_changer.sh -n $n -m $meg_bias -p $i -M $nmeg -P $npeg -d $disp -s $similarity -o $outdir -D $seqDepth -c benchmark_files/benchmark_simul_${param}.txt

		benchmark_imprinting/simulate_genome.sh benchmark_files/benchmark_simul_${param}.txt overwrite #> benchmark_files/log/genome_log.txt
		benchmark_imprinting/simulate_reads.sh benchmark_files/benchmark_simul_${param}.txt #> benchmark_files/log/reads_log.txt
	
		# call imprinting
		run_imprint $param $i
	done
fi

if [ "$falpha" == "true" ]; then
	echo "Changing alpha cutoffs"
	param="alpha"
	array=(0.001 0.005 0.01 0.05 0.15)
	mkdir benchmark_files/$param
	
	printf "Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} similarity: ${similarity}\n" > benchmark_files/$param/wyder_$param.txt
	printf "Anderson method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} similarity: ${similarity}\n" > benchmark_files/$param/anderson_$param.txt
	printf "Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} similarity: ${similarity}\n" > benchmark_files/$param/picard_$param.txt
	printf "Anderson+Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} similarity: ${similarity}\n" > benchmark_files/$param/anderson_picard_combined_$param.txt
	printf "Anderson+Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} similarity: ${similarity}\n" > benchmark_files/$param/anderson_wyder_combined_$param.txt
	printf "Roth method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} similarity: ${similarity}\n" > benchmark_files/$param/roth_$param.txt
	printf "Roth+Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} similarity: ${similarity}\n" > benchmark_files/$param/roth_wyder_combined_$param.txt

	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/picard_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/wyder_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_picard_combined_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_wyder_combined_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/roth_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/roth_wyder_combined_${param}.txt

	benchmark_imprinting/config/read_config_simul.py benchmark_imprinting/config/simulation_config.txt benchmark_files/benchmark_simul_${param}.txt
	benchmark_imprinting/benchmark/config_changer.sh -n $n -M $nmeg -P $npeg -d $disp -m $meg_bias -p $peg_bias -s $similarity -o $outdir -D $seqDepth -c benchmark_files/benchmark_simul_${param}.txt

	benchmark_imprinting/simulate_genome.sh benchmark_files/benchmark_simul_${param}.txt overwrite #> benchmark_files/log/genome_log.txt
	benchmark_imprinting/simulate_reads.sh benchmark_files/benchmark_simul_${param}.txt #> benchmark_files/log/reads_log.txt

	benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_files/benchmark_imprint_${param}.txt

	for i in "${array[@]}"; do
		benchmark_imprinting/benchmark/config_changer.sh -b $mat_bias -B $pat_bias -a $i -o $outdir -l $logfc -C benchmark_files/benchmark_imprint_${param}.txt
		run_imprint $param $i
	done
fi

if [ "$fsim" == "true" ]; then
	echo "Changing %similarity cutoffs"
	param="sim_score"
	array=(80 83 85 87 90 93 95)
	mkdir benchmark_files/$param
	
	printf "Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha}\n" > benchmark_files/$param/wyder_$param.txt
	printf "Anderson method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha}\n" > benchmark_files/$param/anderson_$param.txt
	printf "Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha}\n" > benchmark_files/$param/picard_$param.txt
	printf "Anderson+Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha}\n" > benchmark_files/$param/anderson_picard_combined_$param.txt
	printf "Anderson+Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha}\n" > benchmark_files/$param/anderson_wyder_combined_$param.txt
	printf "Roth method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha}\n" > benchmark_files/$param/roth_$param.txt
	printf "Roth+Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha}\n" > benchmark_files/$param/roth_wyder_combined_$param.txt

	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/picard_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/wyder_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_picard_combined_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_wyder_combined_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/roth_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/roth_wyder_combined_${param}.txt
	
	old_outdir=$outdir
	
	for i in "${array[@]}"; do
		
		outdir=${old_outdir}_${i}
		
		benchmark_imprinting/config/read_config_simul.py benchmark_imprinting/config/simulation_config.txt benchmark_files/benchmark_simul_${param}.txt
		benchmark_imprinting/benchmark/config_changer.sh -n $n -M $nmeg -P $npeg -m $meg_bias -p $peg_bias -d $disp -s $i -o $outdir -D $seqDepth -c benchmark_files/benchmark_simul_${param}.txt

		benchmark_imprinting/simulate_genome.sh benchmark_files/benchmark_simul_${param}.txt overwrite #> benchmark_files/log/genome_log.txt
		benchmark_imprinting/simulate_reads.sh benchmark_files/benchmark_simul_${param}.txt #> benchmark_files/log/reads_log.txt
		
		benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_files/benchmark_imprint_${param}.txt
		benchmark_imprinting/benchmark/config_changer.sh -b $mat_bias -B $pat_bias -a $alpha -o $outdir -l $logfc -C benchmark_files/benchmark_imprint_${param}.txt
	
		run_imprint $param $i
	done
fi

if [ "$fdisp" == "true" ]; then
	echo "Changing dispersion levels"
	param="disp"
	array=("low" "med" "high")
	mkdir benchmark_files/$param
	
	printf "Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/wyder_$param.txt
	printf "Anderson method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/anderson_$param.txt
	printf "Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/picard_$param.txt
	printf "Anderson+Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/anderson_picard_combined_$param.txt
	printf "Anderson+Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/anderson_wyder_combined_$param.txt
	printf "Roth method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/roth_$param.txt
	printf "Roth+Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} similarity: ${similarity}\n" > benchmark_files/$param/roth_wyder_combined_$param.txt

	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/picard_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/wyder_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_picard_combined_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_wyder_combined_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/roth_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/roth_wyder_combined_${param}.txt
	
	for i in "${array[@]}"; do
		
		benchmark_imprinting/config/read_config_simul.py benchmark_imprinting/config/simulation_config.txt benchmark_files/benchmark_simul_${param}.txt
		benchmark_imprinting/benchmark/config_changer.sh -n $n -m $meg_bias -p $peg_bias -M $nmeg -P $npeg -d $i -s $similarity -o $outdir -D $seqDepth -c benchmark_files/benchmark_simul_${param}.txt

		benchmark_imprinting/simulate_genome.sh benchmark_files/benchmark_simul_${param}.txt overwrite #> benchmark_files/log/genome_log.txt
		benchmark_imprinting/simulate_reads.sh benchmark_files/benchmark_simul_${param}.txt #> benchmark_files/log/reads_log.txt
		
		benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_files/benchmark_imprint_${param}.txt
		benchmark_imprinting/benchmark/config_changer.sh -b $mat_bias -B $pat_bias -a $alpha -o $outdir -l $logfc -C benchmark_files/benchmark_imprint_${param}.txt

		# call imprinting
		run_imprint $param $i
	done
fi


if [ "$fseq" == "true" ]; then
	echo "Changing %sequencing depth"
	param="seq_depth"
	array=(0.2 0.25 0.5 1 2 4 5)
	
	mkdir benchmark_files/$param
	
	printf "Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} logfc: ${logfc} similarity: ${similarity}\n" > benchmark_files/$param/wyder_$param.txt
	printf "Anderson method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} logfc: ${logfc} similarity: ${similarity}\n" > benchmark_files/$param/anderson_$param.txt
	printf "Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} logfc: ${logfc} similarity: ${similarity}\n" > benchmark_files/$param/picard_$param.txt
	printf "Anderson+Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} logfc: ${logfc} similarity: ${similarity}\n" > benchmark_files/$param/anderson_picard_combined_$param.txt
	printf "Anderson+Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} logfc: ${logfc} similarity: ${similarity}\n" > benchmark_files/$param/anderson_wyder_combined_$param.txt
	printf "Roth method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} logfc: ${logfc} similarity: ${similarity}\n" > benchmark_files/$param/roth_$param.txt
	printf "Roth+Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${meg_bias} patbias: ${peg_bias} alpha: ${alpha} logfc: ${logfc} similarity: ${similarity}\n" > benchmark_files/$param/roth_wyder_combined_$param.txt
	
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/picard_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/wyder_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_picard_combined_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/anderson_wyder_combined_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/roth_${param}.txt
	printf "${param}\tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark_files/$param/roth_wyder_combined_${param}.txt
	
	old_outdir=$outdir
	
	for i in "${array[@]}"; do
		
		outdir=${old_outdir}_${i}
		
		benchmark_imprinting/config/read_config_simul.py benchmark_imprinting/config/simulation_config.txt benchmark_files/benchmark_simul_${param}.txt
		benchmark_imprinting/benchmark/config_changer.sh -n $n -M $nmeg -P $npeg -m $meg_bias -p $peg_bias -d $disp -s $similarity -o $outdir -D $i -c benchmark_files/benchmark_simul_${param}.txt

		benchmark_imprinting/simulate_genome.sh benchmark_files/benchmark_simul_${param}.txt overwrite #> benchmark_files/log/genome_log.txt
		benchmark_imprinting/simulate_reads.sh benchmark_files/benchmark_simul_${param}.txt #> benchmark_files/log/reads_log.txt
	
		benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_files/benchmark_imprint_${param}.txt
		benchmark_imprinting/benchmark/config_changer.sh -b $mat_bias -B $pat_bias -a $alpha -o $outdir -l $logfc -C benchmark_files/benchmark_imprint_${param}.txt

		run_imprint $param $i
	done
fi
