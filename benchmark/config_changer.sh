#!/bin/bash

while getopts "m:p:M:P:n:r:s:S:d:a:c:C:o:" opt; do
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
		n)	# total_genes = unbiased + megs + pegs
			total_n="$OPTARG"
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
		C)	# config imprint - only needed if changing alpha / outdir
			config_imprint="$OPTARG"
			;;
		o)	# outdir
			outdir="$OPTARG"
			;;
	esac
done

# if [ ${#megs_bias} != 0 ]; then
# 	sed -i "/megs_bias/ { c \
# 	megs_bias=$megs_bias
# 	}" $config_simul
# fi
# 
# if [ ${#pegs_bias} != 0 ]; then
# 	sed -i "/pegs_bias/ { c \
# 	pegs_bias=$pegs_bias
# 	}" $config_simul
# fi

if [ ${#megs} != 0 ]; then
	sed -i "/nmeg/ { c \
	nmeg=$megs
	}" $config_simul
fi

if [ ${#pegs} != 0 ]; then
	sed -i "/npeg/ { c \
	npeg=$pegs
	}" $config_simul
fi

if [ ${#total_n} != 0 ]; then
	sed -i "/total_n/ { c \
	total_n=$total_n
	}" $config_simul
fi

if [ ${#total_n} != 0 ]; then
	unbiased=$(( $total_n - $megs - $pegs ))
	sed -i "/unbiased/ { c \
	unbiased=$unbiased
	}" $config_simul
fi

if [ ${#rep} != 0 ]; then
	sed -i "/rep/ { c \
	rep=$rep
	}" $config_simul
fi

if [ ${#seed} != 0 ]; then
	sed -i "/seed/ { c \
	seed=$seed
	}" $config_simul
fi

if [ ${#disp} != 0 ]; then
	sed -i "/disp/ { c \
	disp=$disp
	}" $config_simul
fi

if [ ${#score} != 0 ]; then
	sed -i "/sim_score/ { c \
	sim_score=$score
	}" $config_simul
fi

if [ ${#outdir} != 0 ]; then
	if [ ${#config_simul} != 0 ]; then
		sed -i "/outdir=/ { c \
		outdir=$outdir
		}" $config_simul
	fi
	if [ ${#config_imprint} != 0 ]; then
		sed -i "/outdir=/ { c \
		outdir=$outdir
		}" $config_imprint
	fi
fi

if [ ${#alpha} != 0 ]; then
	sed -i "/pval_picard/ { c \
	pval_picard=$alpha
	}" $config_imprint
	sed -i "/pval_anderson/ { c \
	pval_anderson=$alpha
	}" $config_imprint
	sed -i "/pval_wyder/ { c \
	pval_wyder=$alpha
	}" $config_imprint
fi

if [ ${#megs_bias} != 0 ]; then
	sed -i "/mat_cutoff_picard/ { c \
	mat_cutoff_picard=$megs_bias
	}" $config_imprint
	sed -i "/mat_cutoff_anderson/ { c \
	mat_cutoff_anderson=$megs_bias
	}" $config_imprint
fi

if [ ${#pegs_bias} != 0 ]; then
	sed -i "/pat_cutoff_picard/ { c \
	pat_cutoff_picard=$pegs_bias
	}" $config_imprint
	sed -i "/pat_cutoff_anderson/ { c \
	pat_cutoff_anderson=$pegs_bias
	}" $config_imprint
fi
