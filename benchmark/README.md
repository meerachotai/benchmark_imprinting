To reproduce benchmarking results in thesis, use `benchmark_imprinting/benchmark` directory

Steps in `benchmarking_pipeline.sh`:
* edits file containing shell environment variables for a few specific changing parameters (using `config_changer.sh`).
* compares the true MEGs/PEGs list and the data analysis pipeline results to write TP/FP table (using `call_true_positives.R`).
* TP/FP table will be saved into `benchmark_files` directory (will be made if it doesn't already exist).

To make figures using the TP/FP table, use `benchmark_imprinting/figures/*.R` files
