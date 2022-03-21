To reproduce benchmarking results in thesis, use `benchmark_imprinting/benchmark` directory

Steps in `benchmarking_pipeline.sh`:
* edits file containing shell environment variables for a few specific changing parameters (using `config_changer.sh`).
* compares the true MEGs/PEGs list and the data analysis pipeline results to write TP/FP table (using `call_true_positives.R`).
* TP/FP table will be saved into `benchmark_files` directory (will be made if it doesn't already exist).

To make figures using the TP/FP table, use `benchmark_imprinting/figures/*.R` files

On the cluster, these commands will run the different benchmarking analyses:
```
qsub -V -N benchmark_pat -cwd -j y -o benchmark_files/log/patbias.txt -m bae -b y -l h_rt=10:00:00,h_data=32G,highp benchmark_imprinting/benchmark/benchmarking_pipeline.sh -p -o benchmark_pat
qsub -V -N benchmark_mat -cwd -j y -o benchmark_files/log/matbias.txt -m bae -b y -l h_rt=10:00:00,h_data=32G,highp benchmark_imprinting/benchmark/benchmarking_pipeline.sh -m -o benchmark_mat
qsub -V -N benchmark_alpha -cwd -j y -o benchmark_files/log/alpha.txt -m bae -b y -l h_rt=10:00:00,h_data=32G,highp benchmark_imprinting/benchmark/benchmarking_pipeline.sh -a -o benchmark_alpha
qsub -V -N benchmark_sim -cwd -j y -o benchmark_files/log/similarity.txt -m bae -b y -l h_rt=10:00:00,h_data=32G,highp benchmark_imprinting/benchmark/benchmarking_pipeline.sh -s -o benchmark_sim
qsub -V -N benchmark_disp -cwd -j y -o benchmark_files/log/disp.txt -m bae -b y -l h_rt=10:00:00,h_data=32G,highp benchmark_imprinting/benchmark/benchmarking_pipeline.sh -d -o benchmark_disp
qsub -V -N benchmark_seqDepth -cwd -j y -o benchmark_files/log/seqDepth.txt -m bae -b y -l h_rt=10:00:00,h_data=32G,highp benchmark_imprinting/benchmark/benchmarking_pipeline.sh -D -o benchmark_seqDepth
```
