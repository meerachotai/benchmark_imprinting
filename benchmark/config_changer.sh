#!/bin/bash

while getopts "m:p:M:P:r:s:S:d:a:c:C:" opt; do
	case $opt in
		m)	# maternal bias simulated
			megs_bias="$OPTARG"
			;;
		p)	# paternal bias simulated
			pegs_bias="$OPTARG"
			;;
		M)	# number of MEGs
			megs="$OPTARG"
			;;
		P)	# number of PEGs
			pegs="$OPTARG"
			;;
		r)	# number of replicates
			rep="$OPTARG"
			;;
		s)	# % similarity
			score="$OPTARG"
			;;
		S)	# seed for reproducibility
			seed="$OPTARG"
			;;
		d)	# dispersion / noise - low/med/high ONLY
			disp="$OPTARG"
			;;
		a)	# for pipelines: p-value / alpha cutoff
			alpha="$OPTARG"
			;;
		c)	# config simul
			config_simul="$OPTARG"
			;;
		C)	# config imprint - only needed if changing alpha
			config_imprint="$OPTARG"
	esac
done

if [ ${#megs_bias} != 0 ]; then
	echo megs_bias=${megs_bias} >> $config_simul
fi

if [ ${#pegs_bias} != 0 ]; then
	echo pegs_bias=${pegs_bias} >> $config_simul
fi

if [ ${#megs} != 0 ]; then
	echo megs=${megs} >> $config_simul
fi

if [ ${#pegs} != 0 ]; then
	echo pegs=${pegs} >> $config_simul
fi

if [ ${#rep} != 0 ]; then
	echo rep=${rep} >> $config_simul
fi

if [ ${#seed} != 0 ]; then
	echo seed=${seed} >> $config_simul
fi

if [ ${#disp} != 0 ]; then
	echo disp=${disp} >> $config_simul
fi

if [ ${#score} != 0 ]; then
	echo score=${score} >> $config_simul
fi

if [ ${#alpha} != 0 ]; then
	echo pval_picard=${alpha} >> $config_imprint
	echo pval_anderson=${alpha} >> $config_imprint
	echo pval_wyder=${alpha} >> $config_imprint
fi
