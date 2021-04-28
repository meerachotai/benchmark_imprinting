# Benchmarking Plant Imprinting Data Analysis Pipelines
## Step 1: Simulating genomes:

Run file: simulate_genome.sh

Note on running file: for now, it does everything in the current directory so navigate to the right directory first (need to change that later)

Dependencies (that are not on loadenv): 
* gffread (https://github.com/gpertea/gffread, should be on $PATH)
* seqkit (https://bioinf.shenwei.me/seqkit/download/)
* simuG (https://github.com/yjx1217/simuG, enter directory for simuG under option -D)

Others (will update): R, samtools

Helper scripts required: (enter directory for helper scripts under option -d)
* edit_genome.sh
* make_annot.R
* inv_transform_sampling.R 

A quick summary of the options is below:
```
-A strainA (for file and folder names)
-B strainB (for file and folder names)
-r original reference genome (from which strainA genome will be simulated)
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
-2 skip to step 2 (only need to rerun step (simulate strainB) when changing similarity %)
-g if there is an error in gffreads, try this option and rerun

(Default score settings as on BLASTn)
```

Sample command:
```
ref="/u/scratch/m/mchotai/rnaseq_simul/col_simul/Cvi.chr.all.v2.0.fasta"
annot="/u/scratch/m/mchotai/rnaseq_simul/col_simul/Cvi.protein-coding.genes.v2.5.2019-10-09.gff3"
scripts_dir="/u/scratch/m/mchotai/rnaseq_simul/scripts_import"
simuG="${scripts_dir}/simuG"
strainA="cviA"
strainB="cviB"

mkdir simul_trial
cd simul_trial

./simulate_genome.sh -r $ref -a $annot -A $strainA -B $strainB -D $simuG -d $scripts_dir -S 70 -s 2 -m 1 -i 0 -e 2.5 -r 1 -n 50 -p 5 -t 10 -T 20 -W 1 -v 5
```
Other notes: I'm currently having some difficulty with placing gffread on $PATH while the script is running. If that is an issue, comment out lines 181-185 and first run gffread directly as given below:
```
mkdir $strainA
gffread -w ${strainA}/${strainA}_transcripts.fa -g $ref $annot
```
Log files of interest:
* per_chrom/un_edited_chr - a list of chromosomes that weren't added to strainB because simuG had to be killed
* per_chrom/indel_scores.txt - stats on the number of trials run per chromosome, including what the final score was (can adjust window accordingly)

The files produced at the end will be: annotations and FASTA files for both strainA and strainB, which can then be used for step 2, where you simulate reads using these simulated genomes.

Note: I have only tried it with 50 chromosomes (-n 50), and it takes around 10 minutes to run through the strainB simulation, might need to optimize runtime 
