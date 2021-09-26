#!/usr/bin/env bash

# Defaults: -------------

disp="med"
mat=95
pat=25
nmeg=200
npeg=200
n=5000
seed=5
alpha=0.05
similarity=90

# Flag options --------------

falpha=false
fdisp=false
fmat=false
fpat=false
fsim=false

while getopts "mpads" opt; do
	case $opt in
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
	esac
done

mkdir benchmark

################
### function ###
################

run_imprint() {
	param=$1
	i=$2
	
	source benchmark_imprinting/config/shell_env_imprint.txt # carry all the variables out, just for $outdir
	
	benchmark_imprinting/simulate_genome.sh benchmark_imprinting/config/shell_env_simul.txt overwrite > benchmark/genome_log.txt
	benchmark_imprinting/simulate_reads.sh benchmark_imprinting/config/shell_env_simul.txt > benchmark/reads_log.txt
	
	benchmark_imprinting/anderson_imprinting.sh benchmark_imprinting/config/shell_env_imprint.txt > benchmark/anderson_log.txt
	
	# for simulation, remove 'strainA'
	awk '{ gsub("strainA", "", $1); print }' $outdir/anderson_MEGs.txt > tmp && mv tmp $outdir/anderson_MEGs.txt
	awk '{ gsub("strainA", "", $1); print }' $outdir/anderson_PEGs.txt > tmp && mv tmp $outdir/anderson_PEGs.txt

	benchmark_imprinting/picard_imprinting.sh benchmark_imprinting/config/shell_env_imprint.txt > benchmark/picard_log.txt
	benchmark_imprinting/wyder_imprinting.sh benchmark_imprinting/config/shell_env_imprint.txt > benchmark/wyder_log.txt

	# output true positive stats
	benchmark_imprinting/benchmark/call_true_positives $outdir > true_pos.txt
	
	printf "tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark/wyder_$param_${i}.txt
		
	meg=$(grep 'total\ Wyder\ maternally-biased:' true_pos.txt | awk '{print $4}')
	true_mat=$(grep 'true\ Wyder\ maternally-biased:' true_pos.txt | awk '{print $4}')
	peg=$(grep 'total\ Wyder\ paternally-biased:'  true_pos.txt | awk '{print $4}')
	true_pat=$(grep 'true\ Wyder\ paternally-biased:' true_pos.txt | awk '{print $4}')

	printf "${totexp}\t$meg\t$peg\t$true_mat\t$true_pat" >> benchmark/wyder_$param_${i}.txt
	
	printf "tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark/anderson_$param_${i}.txt
		
	meg=$(grep 'total\ Anderson\ maternally-biased:' true_pos.txt | awk '{print $4}')
	true_mat=$(grep 'true\ Anderson\ maternally-biased:' true_pos.txt | awk '{print $4}')
	peg=$(grep 'total\ Anderson\ paternally-biased:'  true_pos.txt | awk '{print $4}')
	true_pat=$(grep 'true\ Anderson\ paternally-biased:' true_pos.txt | awk '{print $4}')

	printf "${totexp}\t$meg\t$peg\t$true_mat\t$true_pat" >> benchmark/anderson_$param_${i}.txt
	
	printf "tmat\tpat\ttp_mat\ttp_pat\n" >> benchmark/picard_$param_${i}.txt
		
	meg=$(grep 'total\ Picard\ maternally-biased:' true_pos.txt | awk '{print $4}')
	true_mat=$(grep 'true\ Picard\ maternally-biased:' true_pos.txt | awk '{print $4}')
	peg=$(grep 'total\ Picard\ paternally-biased:'  true_pos.txt | awk '{print $4}')
	true_pat=$(grep 'true\ Picard\ paternally-biased:' true_pos.txt | awk '{print $4}')

	printf "${totexp}\t$meg\t$peg\t$true_mat\t$true_pat" >> benchmark/picard_$param_${i}.txt
}

if [ "$fmat" == "true" ]; then
	echo "Changing %maternal for maternally-biased"
	param="mat"
	array=(75 80 85 90 95 100)
	for i in "${array[@]}"; do
	
		# set up config files
		benchmark_imprinting/config/read_config_simul.py benchmark_imprinting/config/simulation_config.txt benchmark_imprinting/config/shell_env_simul.txt
		benchmark_imprinting/benchmark/config_changer.sh -m $i -p $pat -n $n -M $nmeg -P $npeg -d $disp -s $similarity -c benchmark_imprinting/config/shell_env_simul.txt

		benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_imprinting/config/shell_env_imprint.txt
		benchmark_imprinting/benchmark/config_changer.sh -a $alpha -C benchmark_imprinting/config/shell_env_imprint.txt
	
		printf "Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${i} patbias: ${pat} alpha: ${alpha}\n" > benchmark/wyder_$param_${i}.txt
		printf "Anderson method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${i} patbias: ${pat} alpha: ${alpha}\n" > benchmark/anderson_$param_${i}.txt
		printf "Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${i} patbias: ${pat} alpha: ${alpha}\n" > benchmark/picard_$param_${i}.txt

		# call imprinting
		run_imprint $param $i
	done
fi

if [ "$fpat" == "true" ]; then
	echo "Changing %maternal for paternally-biased"
	param="pat"
	array=(0 10 20 30 40 50)
	for i in "${array[@]}"; do
	
		# set up config files
		benchmark_imprinting/config/read_config_simul.py benchmark_imprinting/config/simulation_config.txt benchmark_imprinting/config/shell_env_simul.txt
		config_changer -m $mat -p $i -n $n -M $nmeg -P $npeg -d $disp -s $similarity -c benchmark_imprinting/config/shell_env_simul.txt

		benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_imprinting/config/shell_env_imprint.txt
		config_changer -a $alpha -C benchmark_imprinting/config/shell_env_imprint.txt
	
		printf "Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${mat} patbias: ${i} alpha: ${alpha} similarity: ${similarity}\n" > benchmark/wyder_$param_${i}.txt
		printf "Anderson method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${mat} patbias: ${i} alpha: ${alpha} similarity: ${similarity}\n" > benchmark/anderson_$param_${i}.txt
		printf "Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${mat} patbias: ${i} alpha: ${alpha} similarity: ${similarity}\n" > benchmark/picard_$param_${i}.txt

		# call imprinting
		run_imprint $param $i
	done
fi

if [ "$falpha" == "true" ]; then
	echo "Changing alpha cutoffs"
	param="alpha"
	array=(0.001 0.005 0.01 0.05)
	for i in "${array[@]}"; do
	
		# set up config files
		benchmark_imprinting/config/read_config_simul.py benchmark_imprinting/config/simulation_config.txt benchmark_imprinting/config/shell_env_simul.txt
		config_changer -m $mat -p $pat -n $n -M $nmeg -P $npeg -d $disp -s $similarity -c benchmark_imprinting/config/shell_env_simul.txt

		benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_imprinting/config/shell_env_imprint.txt
		config_changer -a $i -C benchmark_imprinting/config/shell_env_imprint.txt
	
		printf "Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${mat} patbias: ${pat} alpha: ${i} similarity: ${similarity}\n" > benchmark/wyder_$param_${i}.txt
		printf "Anderson method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${mat} patbias: ${pat} alpha: ${i} similarity: ${similarity}\n" > benchmark/anderson_$param_${i}.txt
		printf "Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${mat} patbias: ${pat} alpha: ${i} similarity: ${similarity}\n" > benchmark/picard_$param_${i}.txt

		# call imprinting
		run_imprint $param $i
	done
fi

if [ "$fsim" == "true" ]; then
	echo "Changing %similarity cutoffs"
	param="similarity"
	array=(80 83 85 87 90 93 95)
	for i in "${array[@]}"; do
	
		# set up config files
		benchmark_imprinting/config/read_config_simul.py benchmark_imprinting/config/simulation_config.txt benchmark_imprinting/config/shell_env_simul.txt
		config_changer -m $mat -p $pat -n $n -M $nmeg -P $npeg -d $disp -s $i -c benchmark_imprinting/config/shell_env_simul.txt

		benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_imprinting/config/shell_env_imprint.txt
		config_changer -a $alpha -C benchmark_imprinting/config/shell_env_imprint.txt
	
		printf "Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${mat} patbias: ${pat} alpha: ${alpha} similarity: ${i}\n" > benchmark/wyder_$param_${i}.txt
		printf "Anderson method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${mat} patbias: ${pat} alpha: ${alpha} similarity: ${i}\n" > benchmark/anderson_$param_${i}.txt
		printf "Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${disp} matbias: ${mat} patbias: ${pat} alpha: ${alpha} similarity: ${i}\n" > benchmark/picard_$param_${i}.txt

		# call imprinting
		run_imprint $param $i
	done

if [ "$fdisp" == "true" ]; then
	echo "Changing dispersion levels"
	param="disp"
	array=("low" "med" "high")
	for i in "${array[@]}"; do
	
		# set up config files
		benchmark_imprinting/config/read_config_simul.py benchmark_imprinting/config/simulation_config.txt benchmark_imprinting/config/shell_env_simul.txt
		config_changer -m $mat -p $pat -n $n -M $nmeg -P $npeg -d $i -s $similarity -c benchmark_imprinting/config/shell_env_simul.txt

		benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_imprinting/config/shell_env_imprint.txt
		config_changer -a $alpha -C benchmark_imprinting/config/shell_env_imprint.txt
	
		printf "Wyder method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${i} matbias: ${mat} patbias: ${pat} alpha: ${alpha} similarity: ${similarity}\n" > benchmark/wyder_$param_${i}.txt
		printf "Anderson method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${i} matbias: ${mat} patbias: ${pat} alpha: ${alpha} similarity: ${similarity}\n" > benchmark/anderson_$param_${i}.txt
		printf "Picard method: n: ${n} nmegs: ${nmeg} npegs: ${npeg} disp: ${i} matbias: ${mat} patbias: ${pat} alpha: ${alpha} similarity: ${similarity}\n" > benchmark/picard_$param_${i}.txt

		# call imprinting
		run_imprint $param $i
	done
fi


