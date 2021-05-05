# Benchmarking Imprinting Data Analysis Pipelines
## Step 1: Simulating genomes

### Run file: simulate_genome.sh

I have only tried it with 50 chromosomes (-n 50), and it takes around 10 minutes to run through the strainB simulation (might need to optimize runtime) 

### Dependencies: 
* gffread (https://github.com/gpertea/gffread, should be on $PATH)
* seqkit (https://bioinf.shenwei.me/seqkit/download/)
* simuG (https://github.com/yjx1217/simuG, enter directory for simuG under option -D)
* Others: R (argparse, tidyverse, Biostrings), samtools

### Helper scripts required: 
(enter directory for helper scripts under option -d)
* edit_genome.sh
* make_annot.R
* inv_transform_sampling.R 
* scoring.R

### Similarity Scoring

A rejection-sampling approach is used for similarity scoring (edit_genome.sh). Given a similarity score, and a SNP-Indel ratio, the number of SNPs and indels that need to be simulated is estimated (scoring.R, -P), which is the input for simuG. The vcf files from simuG are used to calculate the actual score (-V), and compared with desired score. This is repeated until it is within the window of the desired score. This is done per-chromosome basis, which might be slowing down the process. Currently working on modifying this process to a whole-genome basis and pruning the genome as a chromosome achieves the desired score. simuG is intended to work on whole-genome, so this might decrease runtime.

### Options:
```
-o outdirectory (all output will be stored here - HAS to be relative to current working directory)
-A strainA (for file and folder names)
-B strainB (for file and folder names)
-x name for strainA FASTA file
-y name for strainB FASTA file
-X name for strainA annotation file
-Y name for strainB annotation file
-r original reference genome (from which strainA genome will be simulated)
-a original annotation (from which strainA genome will be simulated)
-D simuG directory
-d helper scripts directory
-S similarity score (in %)
-s SNP score (default: 2)
-i Indel score (default: 0)
-e Gap extension score (default: 2.5)
-m Match score (default: 1)
-r SNP:Indel ratio (default: 1)
-n number of 'chromosomes' to be simulated
-v seed
-p time to pause between two simuG runs
-t number of simuG trials to run for each chromosome before killing
-T number of trials of rejection sampling
-W window to allow for rejection sampling (in %)
-2 skip to step 2 (only need to rerun step 2 (simulate strainB) when changing similarity %)
-g if there is an error in gffreads, try this option and rerun

(Default score settings as on BLASTn)
```

Sample command:
```
ref="/path/to/reference.fasta"
annot="/path/to/reference_annotation.gff3"
scripts_dir="/path/to/helper/scripts"
simuG="/path/to/simuG"
simulate_genome.sh -r $ref -a $annot -A strainA -B strainB -x strainA_genome.fa -y strainB_genome.fa -X strainA_annot.gff3 -Y strainB_annot.gff3 -D $simuG -d $scripts_dir -S 70 -s 2 -m 1 -i 0 -e 2.5 -r 1 -n 50 -p 5 -t 10 -T 20 -W 1 -v 5 -o outdir
```
I'm currently having some difficulty with placing gffread on $PATH while the script is running. If that is an issue, comment out lines 181-185 and first run gffread directly as given below:
```
gffread -w strainA/strainA_transcripts.fa -g $ref $annot
```
### Output:

Annotations and FASTA files for both strainA and strainB

Log files of interest:
* $outdir/per_chrom/un_edited_chr.txt - a list of chromosomes that weren't added to strainB because simuG had to be killed
* $outdir/per_chrom/scores_log.txt - stats on the number of trials run per chromosome, including what the final score was (can adjust window accordingly)

## Step 2: Simulating reads

### Run file: simulate_reads.sh

### Dependencies:
* R (argparse, tidyverse, Biostrings)

### Helper scripts required: 
(enter directory for helper scripts under option -d)
* simulated_read-counts.R
* reads_simul.R

### Options:
```
-A strainA (for file and folder names)
-B strainB (for file and folder names)
-x reference strainA FASTA file for read simulation
-y reference strainB FASTA file for read simulation
-d helper scripts directory
-m number of MEGs to be simulated
-p number of PEGs to be simulated
-u number of unbiased genes to be simulated
-M %maternal bias for MEGs
-P %maternal bias for PEGs
-r read length for FASTQ files
-R number of replicates to be simulated (needed for DESeq2 specifically)
-s seed
-o outdirectory (HAS to be relative to current directory, will store output files here)
```

Sample command:
 ```
scripts_dir="/path/to/helper_scripts/scripts_import"
refA="/path/to/refA.fa"
refB="/path/to/refB.fa"

$scripts_dir/simulate_reads_opt.sh -A strainA -B strainB -x $refA -y $refB -d $scripts_dir -s 5 -u 30 -m 10 -p 10 -r 50 -R 3 -M 95 -P 25 -o outdir
```

### Output:

FASTQ (.fq) files that match counts simulated.

Other relevant files: 
* $outdir/reads_simul/simul_counts+id_A.txt and $outdir/reads_simul/simul_counts+id_B.txt - a summary of 'chromosome' ids alongside read counts
* $outdir/counts_simul_megs.txt and $outdir/counts_simul_pegs.txt - true MEG and PEG lists to use in verifying imprinting calls


