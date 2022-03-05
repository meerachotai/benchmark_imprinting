# param=counts
# source benchmark_files/benchmark_imprint_${param}.txt

source $1

logfc=${logfc_wyder}
pval=${pval_wyder}
outprefix="wyder"
counts_dir=${counts_dir_wyder}

map=$outdir/map
count=$rep
outprefix=anderson_wyder


$scripts_dir/call_imprinting_wyder.R -C -c "${map}/anderson_" -r $count -p $pval -l $logfc ${outdir}/$outprefix
