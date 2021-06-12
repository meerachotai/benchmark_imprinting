# Benchmarking Imprinting Data Analysis Pipelines
## Step 1: Simulating genomes

### Run file: `simulate_genome.sh`

### Dependencies: 
* gffread (https://github.com/gpertea/gffread, should be on $PATH)
* seqkit (https://bioinf.shenwei.me/seqkit/download/)
* simuG (https://github.com/yjx1217/simuG, enter directory for simuG under option -D)
* samtools
* Others: R (argparse, tidyverse, Biostrings)

### Helper scripts required: 
(enter directory for helper scripts under option -d)
* `edit_genome.sh`
* `make_annot.R`
* `inv_transform_sampling.R`
* `scoring.R`

### Similarity Scoring

A rejection-sampling approach is used for similarity scoring (`edit_genome.sh`). Given a similarity score, and a SNP-Indel ratio, the number of SNPs and indels that need to be simulated is estimated (`scoring.R`, -P), which is the input for simuG. The vcf files from simuG are used to calculate the actual score (-V), and compared with desired score. This is repeated until it is within the window of the desired score. This is done per-chromosome basis in the main branch given below.

**Alternate branch (`helper_scripts/edit_genome_v2`):** Operates on a whole-genome basis and prunes the original genome as a chromosome achieves the desired score. This increased runtime relative to the original per-chromosome method (see metrics below). Increasing the window might decrease runtime, possibly, but it hasn't been tried out yet.

Added dependencies for alternate scripts: R - vcfR library

Sample command (after conducting step 1 from below):
```
edit_genome_v2.sh -A $strainA -a $refA -B $strainB -b refB -D $simuG_dir -d $scripts_dir -S 90 -n 50 -v 5 -p 3 -t 10 -W 3 -o $outdir
```

Some one-time metrics on comparing the two approaches:
* per-chromosome basis ran for 40 minutes for 50 chromosomes
* whole-genome basis ran for 3 hours for 39 chromosomes

### Options:
```
-o outdirectory (all output will be stored here - HAS to be relative to current working directory)
-A strainA (for file and folder names)
-B strainB (for file and folder names)
-x outprefix for strainA FASTA file
-y outprefix for strainB FASTA file
-X outprefix for strainA annotation file
-Y outprefix for strainB annotation file
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
simulate_genome.sh -r $ref -a $annot -A strainA -B strainB -x strainA_genome -y strainB_genome -X strainA_annot -Y strainB_annot -D $simuG_dir -d $scripts_dir -S 70 -s 2 -m 1 -i 0 -e 2.5 -r 1 -n 50 -p 5 -t 10 -T 20 -W 1 -v 5 -o outdir
```
I'm currently having some difficulty with placing gffread on $PATH while the script is running. If that is an issue, comment out lines 205-209 and first run gffread directly as given below:
```
gffread -w outdir/strainA/strainA_transcripts.fa -g $ref $annot
```
### Output:

Annotations and FASTA files for both strainA and strainB (in outdir)

Log files of interest:
* `$outdir/per_chrom/un_edited_chr.txt` - a list of chromosomes that weren't added to strainB because simuG had to be killed
* `$outdir/per_chrom/scores_log.txt` - stats on the number of trials run per chromosome, including what the final score was (can adjust window accordingly)

## Step 2: Simulating reads

### Run file: `simulate_reads.sh`

### Dependencies:
* R (argparse, tidyverse, Biostrings)

### Helper scripts required: 
(enter directory for helper scripts under option -d)
* `simulated_read-counts.R`
* `reads_simul.R`

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
simulate_reads.sh -A strainA -B strainB -x $refA -y $refB -d $scripts_dir -s 5 -u 30 -m 10 -p 10 -r 50 -R 3 -M 95 -P 25 -o outdir
```

### Output:

FASTQ (.fq) files that match counts simulated.

Other relevant files: 
* `$outdir/reads_simul/simul_counts+id_A.txt` and `$outdir/reads_simul/simul_counts+id_B.txt` - a summary of 'chromosome' ids alongside read counts
* `$outdir/counts_simul_megs.txt` and `$outdir/counts_simul_pegs.txt` - true MEG and PEG lists to use in verifying imprinting calls

## Step 3: Calling Imprinting

Steps 1 and 2 can be skipped in favor of providing real-data files for mapping and calling imprinting Step 3 onwards. However, there are some file name conventions that **must** be followed in order for this to work - refer to the next section below. 

Depending on the method being used, some additional conventions are also required for input files, so refer to individual Method sections for more information.

### Required conventions:

#### FASTQ files
For this particular step, a specific FASTQ file name is required in order to map the reads correctly.

* The file name MUST end with the suffix **`cross_replicate.fq`**
* The address and the prefix of the file names should be placed under the option -f

**Example:** 

The simulation steps above have the following file for the first (1) replicate, for the replicate cross AxB with the prefix being `$strainA_$strainB_`. Adding the file's location to the prefix, the file is called: `$outdir/reads_simul/$strainA_$strainB_AxB_1.fq`.

Following the above convention, under -f option place: -f `$outdir/reads_simul/$strainA_$strainB_`

If you followed with the steps 1 and 2, there's no need to use -f at all, it is done for you by default.

## Method A: Anderson et al. 

### Run file: `anderson_mapping.sh`

Adapted from scripts on: https://github.com/SNAnderson/Imprinting2020

### Dependencies:
* hisat2
* HTSeqcount
* R (DESeq2, argparse)

### Helper scripts required:
(enter directory for helper scripts under option -d)
* `call_imprinting_anderson.R`

### Required conventions:

#### FASTQ files
For this particular step, a specific FASTQ file name is required in order to map the reads correctly.

* The file name MUST end with the suffix **`cross_replicate.fq`**
* The address and the prefix of the file names should be placed under the option -f

**Example:** 

The simulation steps above have the following file for the first (1) replicate, for the replicate cross AxB with the prefix being `$strainA_$strainB_`. Adding the file's location to the prefix, the file is called: `$outdir/reads_simul/$strainA_$strainB_AxB_1.fq`.

Following the above convention, under -f option place: -f `$outdir/reads_simul/$strainA_$strainB_`

If you followed with the steps 1 and 2, there's no need to use -f at all, it is done for you by default.

#### FASTA and annotation files

For the annotation files, depending on your file's convention, the 'attributes' column will have a different label (ID, Parent etc.). HTseqcount needs to know what this label is under its option -i. This script also has an option -i for you to provide this, with the default being set as ID. For more information on the conventions, refer: http://gmod.org/wiki/GFF3#GFF3_Format

If you followed with steps 1 and 2, enter the annotation with the `_anderson.gff3` suffix under the -X and -Y options. Additionally, there is no need to use the -i option, it uses the default.

Since the Anderson method concatenates reference and annotation files before read mapping, you are required to make sure the chromosome names for different strains are different. For this, a function is set up within the script that renames the chromosome names in the files according to strain names, but note that this does NOT apply to all reference and annotation file conventions. 

The renaming function is used automatically, but if you want to skip it, use the option -e. If you followed with steps 1 and 2, do not use -e, and it will work with your files reference and annotation files automatically if you provide the correct location for it under the -x -y -X and -Y options.

Below is the function does with some comments and preview of what it does, which might help if you're trying to edit it beforehand.

```
# strain is the strain name for A/B strain (ex: cviA)
# ref is the reference FASTA file that needs to be edited
# annot is the annotation .gff3 file that needs to be edited
 
# for all lines the begin with >, replace > with >$strain_
sed "s/^>/>${strain}_/g" $ref > ${strain}_genome.fa 

# at the start of each line, it appends with prefix $strain
sed "s/^/${strain}_/" $annot > map/${strain}_annot.gff3 # start of line	

# since the previous command also applied on the 'gff-version 3' header in the annotation file, replace it	
head="##gff-version 3"
sed -i "1s/.*/$head/" map/${strain}_annot.gff3
```

**FASTA:** 

Before:
```
>ATCVI-1G19970.1 gene=ATCVI-1G19970 CDS=1-1998
ATGTTCATGTCTGTGTTTGAAATTGCATTCTGCCAAAACCAAACTCCTGAGCCAGAATCA
ACACAAGTCTTGCAGCTACATTACCAAGATCGATATGTGGTAGTGGAGAATAGATTTTTA
CAATTAACTTTATCAAATCCTGAGGGTTTCGTCACCGGAATCCAGTATAATGGTATCGAC 
```
After:
```
>cviA_ATCVI-1G19970.1 gene=ATCVI-1G19970 CDS=1-1998
ATGTTCATGTCTGTGTTTGAAATTGCATTCTGCCAAAACCAAACTCCTGAGCCAGAATCA
ACACAAGTCTTGCAGCTACATTACCAAGATCGATATGTGGTAGTGGAGAATAGATTTTTA
CAATTAACTTTATCAAATCCTGAGGGTTTCGTCACCGGAATCCAGTATAATGGTATCGAC
```
**Annotation:**

Before:
```
##gff-version 3
ATCVI-1G19970.1	.	exon	1	1961	.	+	.	ID=ATCVI-1G19970.1A
ATCVI-1G33940.1	.	exon	1	1544	.	+	.	ID=ATCVI-1G33940.1A
ATCVI-1G38040.1	.	exon	1	401	.	+	.	ID=ATCVI-1G38040.1A
ATCVI-1G42830.1	.	exon	1	1141	.	+	.	ID=ATCVI-1G42830.1A
```
After:
```
##gff-version 3
cviA_ATCVI-1G19970.1	.	exon	1	1961	.	+	.	ID=ATCVI-1G19970.1A
cviA_ATCVI-1G33940.1	.	exon	1	1544	.	+	.	ID=ATCVI-1G33940.1A
cviA_ATCVI-1G38040.1	.	exon	1	401	.	+	.	ID=ATCVI-1G38040.1A
cviA_ATCVI-1G42830.1	.	exon	1	1141	.	+	.	ID=ATCVI-1G42830.1A
```
#### Gene Key

If steps 1 and 2 were not followed, a gene key is **required** for calling imprinting. It must be a tab-delimited file .txt with a header. It should consist of a list of syntelogs between strainA and strainB. Make sure that the syntelog names match with the 'attribute' column of the annotation files. It must have strainA as first column, and strainB as the second column. Enter gene key filename under option -g.

**Example:**
```
A B
ATCVI-1G19970.1A	ATCVI-1G19970.1B
ATCVI-1G33940.1A	ATCVI-1G33940.1B
ATCVI-1G38040.1A	ATCVI-1G38040.1B
ATCVI-1G42830.1A	ATCVI-1G42830.1B
ATCVI-1G64390.1A	ATCVI-1G64390.1B
ATCVI-1G64620.1A	ATCVI-1G64620.1B
ATCVI-1G67300.1A	ATCVI-1G67300.1B
ATCVI-1G75310.1A	ATCVI-1G75310.1B
ATCVI-1G81680.1A	ATCVI-1G81680.1B
```
#### Counts Files

If you already have the counts files and are planning to use the helper script `call_imprinting_anderson.R` directly, use one of the two options below:
* name files exactly with the suffix **`cross_replicate.txt`**, provide the prefix under the option -c AND use option -C to indicate that the files need to be concatenated
* merge the files with the first 1:replicates columns for the AxB counts, and replicates+1:replicates\*2 columns for the BxA counts, and provide the one filename under -c AND DO NOT use the option -C to indicate that the file is already concatenated. An example for a concatenated file with 3 replicates is given below:

```
                 AxB_1 AxB_2 AxB_3 BxA_1 BxA_2 BxA_3
ATCVI-1G19970.1A   228    74    38     6    38     7
ATCVI-1G19970.1B    35     7    14    90   202    89
ATCVI-1G33940.1A    48   299   110    65    22    50
ATCVI-1G33940.1B   135    46    53   256    69   126
ATCVI-1G38040.1A    55   156   109    91    60    17
ATCVI-1G38040.1B    99    63    11   189    85   159
```

Note that for Anderson/DESeq2 approach, replicates are required.

### Options:
```
-o outdirectory (all output will be stored here - HAS to be relative to current working directory)
-A strainA (for file and folder names)
-B strainB (for file and folder names)
-x name for strainA FASTA file
-y name for strainB FASTA file
-X name for strainA annotation file
-Y name for strainB annotation file
-d helper scripts directory
-r number of replicates available
-i annotation gff3 attribute name (column 9) (default: ID)
-f represents the start of the FASTQ reads file name
-e boolean, use if you want to skip editing the files
-p p-value / alpha cutoff for calling imprinting (DESeq2's padj) (default: 0.05)
-l log2fc cutoff for calling imprinting (DESeq2's log2foldchange) (default: 1, recommended by authors)
-M maternal expression cutoff for calling imprinting (default: 0.8)
-P paternal expression cutoff for calling imprinting
-g gene key with syntelogs between strainA and strainB
-a attribute substring that is UNIQUE to strainA
-b attribute substring that is UNIQUE to strainB
-O outprefix for imprinting files
```

Sample command:
```
anderson_mapping.sh -A $strainA -B $strainB -x $refA -y $refB -X $annotA -Y $annotB -o $outdir -e -i ID -r 3 -d $scripts_dir -f $fastq_dir -O $outdir/outprefix -a .1A -b .1B -M 0.8 -P 0.5 -l 1 -p 0.05 -g gene_key.txt
```

### Output:

* `outprefix_anderson_MEGs.txt` and `outprefix_anderson_PEGs.txt` - imprinted MEGs/PEGs lists with both syntelog names given
* `outprefix_anderson_stats.txt` - DESeq2 stats and imprinting status for all syntelogs

## Method B: Picard and Gehring 

### Run file: `picard_mapping.sh`

Wrapper around the suite of scripts on: https://github.com/clp90/imprinting_analysis. For more control over parameters, use this suite directly.

### Dependencies:
* Picard and Gehring's Imprinting Analysis suite (https://github.com/clp90/imprinting_analysis, enter directory under -p). Must also include all its dependencies, as outlined on their README file.

### Options:
```
-o outdirectory (all output will be stored here - HAS to be relative to current working directory)
-p Picard and Gehring's Imprinting analysis suite directory
-A strainA (for file and folder names)
-B strainB (for file and folder names)
-a reference/strainA annotation file
-g reference/strainA genome FASTA file
-M MEGs cutoff for calling imprinting (default: 95)
-P PEGs cutoff for calling imprinting (default: 25)
-s SNPs file, in required format
-r number of replicates (default: 3)
-f represents the start of the FASTQ reads file name
```
### Required Conventions:

#### For Picard and Gehring's suite
Some conventions that need to be followed for the program to run correctly:
* Annotations file must have transcript_id and gene_id under attributes column (for STAR and HTSeqcount, respectively). 
* Annotation files should not contain the ##gff-version 3 header for STAR
* SNPs file must be in the required format. 
* Names for chromosomes in genome FASTA files/read files/annotation files should not contain a '.'

If you followed with steps 1 and 2, adjustments will be made automatically:
* Use the annotation file with the suffix `_picard.gff3` under the -a option.
* Leave the -s option blank, a SNP file will be made from the vcf files under outdir/per_chrom. 

#### FASTQ files
For this particular step, a specific FASTQ file name is required in order to map the reads correctly.

* The file name MUST end with the suffix **`cross_replicate.fq`**
* The address and the prefix of the file names should be placed under the option -f

**Example:** 

The simulation steps above have the following file for the first (1) replicate, for the replicate cross AxB with the prefix being `$strainA_$strainB_`. Adding the file's location to the prefix, the file is called: `$outdir/reads_simul/$strainA_$strainB_AxB_1.fq`.

Following the above convention, under -f option place: -f `$outdir/reads_simul/$strainA_$strainB_`

If you followed with the steps 1 and 2, there's no need to use -f at all, it is done for you by default.

For now, FASTQ files from reciprocal crosses with same replicate number attached to them are considered paired. Working on adding another option to consider them separately (calling imprinting using combinations of all AxB and BxA replicates).

Sample command:

```
picard_mapping.sh -A strainA -B strainB -g $refA -a $annotA -r 3 -M 95 -P 25 -p $picard -o $outdir
```
### Output:

* `outdir/picard_map/rep_x_imprinting/imprinting/rep_x_imprinting_filtered_MEGs.txt` and `outdir/picard_map/rep_x_imprinting/imprinting/rep_x_imprinting_filtered_PEGs.txt` where 'x' represents replicate number.

