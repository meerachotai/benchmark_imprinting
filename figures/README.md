To reproduce figures in thesis, use `benchmark_imprinting/figures` directory

Requires TP/FP tables and some specific outdirectories written when running `benchmarking_pipeline.sh`. 

Scripts and the figures that they produce:
* Fig. 3: `anderson_false_neg.R` 
* Fig. 4 and 7: `anderson_picard_counts.R`
* Fig. 5-6 and 8-9: `benchmarking_plots.R`
* Fig. 10: `benchmarking_seq_depth.R`

More details on scripts and inputs follow:

**Fig. 3** requires:
* `anderson_stats.txt`
* `true_*.txt` 

from one of the outdirectories (for default settings, I used `benchmark_sim_90`)

**Fig 4** requires:
* `benchmark_sim_90/counts` directory
* `benchmark_sim_90/simul_counts+id_*.txt`

For **Fig. 7**, use `benchmark_sim_95` instead

**Fig. 5-6 and 8-9** require:
* `benchmark_files` directory

**Fig. 10** requires:
* `benchmark_files/seq_depth` directory
* `benchmark_seqDepth_*` directories
* `benchmark_files/seqDepth_counts` directory (for instructions on how to get the contents of this directory from `benchmark_seqDepth_*`, refer to code-block below)

```
array=(0.2 0.25 0.5 1 2 4 5)
rep=3

mkdir benchmark_files/seqDepth_counts

# copy Picard-generated counts
for i in "${array[@]}"; do
  for j in $(seq 1 1 $rep); do
    cp benchmark_seqDepth_${i}/picard/picard_map/rep_${j}_${j}_imprinting/counts_per_gene/rep_${j}_${j}_strainAxstrainB_counts_merged.txt seqDepth_counts/AxB_${i}_${j}_counts.txt
  done
done

# get lengths of each 'chromosome' in transcriptome
for i in "${array[@]}"; do
	cut -f1-2 benchmark_seqDepth_${i}/ref_A.fa.fai > seqDepth_counts/length_id_${i}.txt
done
```

