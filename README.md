# Benchmarking Imprinting Data Analysis Pipelines
## Step 1: Simulating genomes

### Dependencies: 
* gffread (https://github.com/gpertea/gffread)
* samtools
* Python (> 3; argparse, scipy, random, numpy)
* R (argparse, tidyverse, Biostrings)

### Helper scripts required: 
(enter directory for helper scripts under option -d)
* `edit_genome.py`
* `make_annot.R`

### Run:
Fill in `benchmark_imprinting/config/simulation_config.txt` as per your requirements.

Run the following scripts: 
* `benchmark_imprinting/config/read_config_simul.py benchmark_imprinting/config/simulation_config.txt benchmark_imprinting/config/shell_env_simul.txt`
* `benchmark_imprinting/simulate_genome.sh benchmark_imprinting/config/shell_env_simul.txt`

### Output:

* annotations and FASTA files for both strainA and strainB in the out-directory.
* `$outdir/all_genes.txt` - list of genes selected from gffread-transcripts generated
* `$strainB/$strainB_ref2sim.txt` - reference vs. simulated genome for SNPs and indels. 
* `$strainB/$strainB_similarity.txt` - similarity % of each individual chromosome

## Step 2: Simulating reads

### Dependencies:
* R (argparse, tidyverse, Biostrings)

### Helper scripts required: 
(enter directory for helper scripts under option -d)
* `simulate_counts.R`
* `reads_simul.R`

### Run:
`benchmark_imprinting/simulate_reads.sh benchmark_imprinting/config/shell_env_simul.txt`

### Output:

FASTQ (.fq) files that match counts simulated.

Other relevant files: 
* `$outdir/reads_simul/simul_counts+id_A.txt` and `$outdir/reads_simul/simul_counts+id_B.txt` - a summary of gene ids alongside read counts
* `$outdir/true_MEGs.txt` and `$outdir/true_PEGs.txt` - MEGs and PEGs lists to use in verifying imprinting calls

## Step 3: Calling Imprinting

Steps 1 and 2 can be skipped in favor of providing real-data files for mapping and calling imprinting Step 3 onwards. However, there are some file name conventions that **must** be followed in order for these methods to work - refer to individual Method sections for more information.

## Method A: Anderson et al. 

Adapted from scripts on: https://github.com/SNAnderson/Imprinting2020

### Dependencies:
* hisat2
* HTSeqcount
* R (DESeq2, argparse)

### Helper scripts required:
(enter directory for helper scripts under option -d)
* `get_counts_anderson.R`
* `call_imprinting_anderson.R`

### Required conventions:

#### FASTQ files
For this particular step, a specific FASTQ file name is required in order to map the reads correctly.

* The file name MUST end with the suffix **`cross_replicate.fq`**
* The address and the prefix of the file names should be placed under FASTQ_DIR
* If the reads are paired-ended, the end of the file's name should be `_1.fq` and `_2.fq`

**Example:** 

The simulation steps above have the following file for the first (1) replicate, for the replicate cross AxB with the prefix being `$strainA_$strainB_`. Adding the file's location to the prefix, the file is called: `$outdir/reads_simul/$strainA_$strainB_AxB_1.fq`.

Following the above convention, your config file should read: `FASTQ_DIR = $outdir/reads_simul/$strainA_$strainB_`

Alternately - provide a tab-delimited file with columns in the following order:
1. Filename (must be FASTQ format, must have extension \*.txt, \*.fastq, or \*.fq)
2. Sample's replicate number (should only go up to $rep)
3. 'F' if single-end reads / forward-paired end reads or 'R' if reverse paired-end reads
4. Sample's reciprocal cross - 'AxB' or 'BxA'

**Example:** see `config/single-end_config.txt` or `config/paired-end_config.txt`

#### FASTA and annotation files

For the annotation files, depending on your file's convention, the 'attributes' column will have a different label (ID, Parent etc.). HTSeqcount needs to know what this label is under its option -i. Use the HTSEQ_OPTION_I in the config file to provide this information. For more information on the conventions, refer: http://gmod.org/wiki/GFF3#GFF3_Format.

If you followed with steps 1 and 2, enter the annotation with the `_anderson.gff3` suffix under the Anderson annotation files in the config file and set `HTSEQ_OPTION_I = ID`.

Since the Anderson method concatenates reference and annotation files before read mapping, you are required to make sure the chromosome names for different strains are different. For this, a function is set up within the script that renames the chromosome names in the files according to strain names, but note that this does NOT apply to all reference and annotation file conventions. To use the renaming function, set `RENAME_HEADERS = TRUE`.

#### Gene Key

If steps 1 and 2 were not followed, a gene key is **required** for calling imprinting. It must be a tab-delimited .txt file with a header. It should consist of a list of syntelogs between strainA and strainB. Make sure that the syntelog names match with the 'attribute' column of the annotation files. It must have strainA as first column, and strainB as the second column. Enter gene key filename under `GENE_KEY`. Under `UNIQUE_ANNOT_ATTRIBUTE_ID_A/B`, enter the substrings that uniquely belong to the A and B gene names respectively. 

**Example:** gene_key.txt (the header does not necessarily have to be the one given below:)
```
A B
ATCVI-1G19970cviA ATCVI-1G19970cviB
ATCVI-1G33940cviA ATCVI-1G33940cviB
ATCVI-1G38040cviA ATCVI-1G38040cviB
ATCVI-1G42830cviA ATCVI-1G42830cviB
ATCVI-1G64390cviA ATCVI-1G64390cviB
ATCVI-1G64620cviA ATCVI-1G64620cviB
ATCVI-1G67300cviA ATCVI-1G67300cviB
ATCVI-1G75310cviA ATCVI-1G75310cviB
ATCVI-1G81680cviA ATCVI-1G81680cviB
```

For the example above, enter: 
```
GENE_KEY = gene_key.txt
UNIQUE_ANNOT_ATTRIBUTE_ID_A cviA
UNIQUE_ANNOT_ATTRIBUTE_ID_B cviB
```

#### Counts Files

If you already have the counts files and are planning to use the helper script `call_imprinting_anderson.R` directly, use one of the two options below:
* Name files exactly with the suffix **`cross_replicate.txt`**, provide the prefix under the option -c AND use option -C to indicate that the files need to be concatenated. Each file should have two columns: first column with gene names (including genes from strainA and strainB, labelled uniquely, as given under options -a and -b) and second column with counts.

**Example:** `counts_cviA_cviB_AxB_1.txt` would require that your input command includes: `-c counts_cviA_cviB_ -C`

* Merge the files with a header, the gene names as rownames, with alternating AxB and BxA columns. Provide the one filename under -c and do NOT use the option -C to indicate that the file is already concatenated. An example for a concatenated file with 3 replicates is given below:

**Example:** (the header does not necessarily have to be the one given below:)
```
                 AxB_1 BxA_1 AxB_2 BxA_2 AxB_3 BxA_3
ATCVI-1G19970cviA   228   74    38     6    38     7
ATCVI-1G19970cviB    35    7    14    90   202    89
ATCVI-1G33940cviA    48  299   110    65    22    50
ATCVI-1G33940cviB   135   46    53   256    69   126
ATCVI-1G38040cviA    55  156   109    91    60    17
ATCVI-1G38040cviB    99   63    11   189    85   159
```

Note that for Anderson/DESeq2 approach, > 1 replicates are required.

### Run:

Fill in `benchmark_imprinting/config/imprinting_config.txt` as per your requirements.

Run the following scripts: 
* `benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_imprinting/config/shell_env_imprint.txt`
* `benchmark_imprinting/anderson_imprinting.sh benchmark_imprinting/config/shell_env_imprint.txt`

### Output:

* `anderson_MEGs.txt` and `anderson_PEGs.txt` - imprinted MEGs/PEGs lists with both syntelog names given
* `anderson_stats.txt` - DESeq2 stats and imprinting status for all syntelogs

Output from `get_counts_anderson.R` is in a similar format to the output from the Picard pipeline, so that the second half of the other pipelines (Picard, Wyder) can be used to call imprinting. The file names are as below:
* `anderson_cross_rep.txt`files have columns with gene-names, A counts, B counts.
* `anderson_cross_A_rep.txt` and `outprefix_cross_B_rep.txt` have the A and B counts separately as well.

## Method B: Picard and Gehring 

Wrapper around the suite of scripts on: https://github.com/clp90/imprinting_analysis. For more control over parameters, use this suite directly. This wrapper script also has the added functionality of running all reciprocal replicate pairs or all combinations of pairs of reciprocal replicates to call imprinted genes, and then consensus calling to give a final list of imprinted genes.

### Dependencies:
* Picard and Gehring's Imprinting Analysis suite (https://github.com/clp90/imprinting_analysis, enter directory under `PICARD_SCRIPTS_DIRECTORY`). Must also include all its dependencies, as outlined on their README file.

### Helper scripts required:
(enter directory for helper scripts under option -d)
* `find_consensus.py`

Other notes: this wrapper script requires a large memory, and depending on the size of your genome, and the number of replicates available it could take very long to run as well - set up jobs/interactive sessions accordingly.

### Required Conventions:

#### FASTQ files
For this particular step, a specific FASTQ file name is required in order to map the reads correctly.

* The file name MUST end with the suffix **`cross_replicate.fq`**
* The address and the prefix of the file names should be placed under the option -f

**Example:** 

The simulation steps above have the following file for the first (1) replicate, for the replicate cross AxB with the prefix being `$strainA_$strainB_`. Adding the file's location to the prefix, the file is called: `$outdir/reads_simul/$strainA_$strainB_AxB_1.fq`.

Following the above convention, your config file should read: `FASTQ_DIR = $outdir/reads_simul/$strainA_$strainB_`

Alternately - provide a tab-delimited file with columns in the following order:
1. Filename (must be FASTQ format, must have extension \*.txt, \*.fastq, or \*.fq)
2. Sample's replicate number (should only go up to $rep)
3. 'F' if single-end reads / forward-paired end reads or 'R' if reverse paired-end reads
4. Sample's reciprocal cross - 'AxB' or 'BxA'

**Example:** see `config/single-end_config.txt` or `config/paired-end_config.txt`

#### For Picard and Gehring's suite
Some conventions that need to be followed for the program to run correctly:
* SNPs file must be in the required format. 
* Annotations file must have transcript_id and gene_id under attributes column (for STAR and HTSeqcount, respectively). 
* Annotation files should not contain the ##gff-version 3 header for STAR
* Names for chromosomes in genome FASTA files/read files/annotation files should not contain a '.'

An example is given below:
```
ATCVI-1G19970	.	exon	1	1998	.	+	.	gene_id=ATCVI-1G19970;transcript_id ATCVI-1G19970cviA
ATCVI-1G33940	.	exon	1	1548	.	+	.	gene_id=ATCVI-1G33940;transcript_id ATCVI-1G33940cviA
ATCVI-1G38040	.	exon	1	411	.	+	.	gene_id=ATCVI-1G38040;transcript_id ATCVI-1G38040cviA
ATCVI-1G42830	.	exon	1	1134	.	+	.	gene_id=ATCVI-1G42830;transcript_id ATCVI-1G42830cviA
ATCVI-1G64390	.	exon	1	1524	.	+	.	gene_id=ATCVI-1G64390;transcript_id ATCVI-1G64390cviA
ATCVI-1G64620	.	exon	1	1041	.	+	.	gene_id=ATCVI-1G64620;transcript_id ATCVI-1G64620cviA
```

If you followed with steps 1 and 2:
* Set `ANNOT_A = $annotA_picard.gff3`
* Leave the SNPS_FILE option blank, a SNP file will be made from the vcf files under outdir/per_chrom. 

If you have a different number of AxB and BxA replicates, specifically add these in under AxB_REP and BxA_REP respectively, so that combinations of the two sets will be used to find consensus MEGs and PEGs.

### Run:

Fill in `benchmark_imprinting/config/imprinting_config.txt` as per your requirements.

Run the following scripts: 
* `benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_imprinting/config/shell_env_imprint.txt`
* `benchmark_imprinting/picard_imprinting.sh benchmark_imprinting/config/shell_env_imprint.txt`

### Output:

* `picard_all_MEGs.txt` and `picard_all_PEGs.txt` - all MEGs/PEGs list, along with combination count from which it was called
* `picard_MEGs.txt` and `picard_PEGs.txt` - consensus MEGs/PEGs list
* `${outdir}/picard_map/rep_${i}_${j}_imprinting/counts_per_gene` - directory with counts files generated (where i and j are AxB and BxA replicate numbers respectively)

## Method C: Wyder et al.

Adapted from scripts on: https://github.com/swyder/Reanalysis_plant_imprinting

### Dependencies:
* R (argparse, edgeR)

### Helper scripts required:
(enter directory for helper scripts under option -d)
* `call_imprinting_wyder.R`

### Required Conventions:

#### Counts Files

Running `wyder_imprinting.sh` assumes that you have first run `picard_mapping.sh` and have acquired the `rep_${i}_${j}_imprinting/counts_per_gene/rep_${i}_${j}_${strainA}x${strainB}_counts_merged.txt` files already. By default it assumes that these counts files are under the directory: `$outdir/picard_map.` Provide an alternate prefix for these counts files under `COUNTS_DIRECTORY`. Maintain the same settings for paired replicates. 

You can also instead run the mapping and counting steps recommended by Wyder et al. method at https://github.com/swyder/Reanalysis_plant_imprinting, and then directly use the call_imprinting_wyder.R script to run edgeR to call imprinted genes. Any alternate method for mapping and counting also works. However, this will mean that the counts files will need to be in one of two formats given below:
* Name files exactly with the suffix **`cross_replicate.txt`**, provide the prefix under the option -c AND use option -C to indicate that the files need to be concatenated. Each file should have three columns: the gene name, strainA counts and strainB counts.

**Example:** `wyder_AxB_1.txt` would require that your input command includes: `-c wyder_ -C`

* Merge the files with a header, gene names as row names, the first 1:(2\*replicates) columns for the AxB counts, and (2\*replicates)+1:(replicates\*4) columns for the BxA counts (with strainA counts first, and strainB counts second). Provide the one filename under -c and do NOT use the option -C to indicate that the file is already concatenated. An example for a concatenated file with 1 replicate is given below:

**Example:** (the header does not necessarily have to be the one given below:)
```
                 AxB_1_A AxB_1_B BxA_1_A BxA_1_B 
ATCVI-1G19970     228   74    38     6    
ATCVI-1G33940      48  299   110    65  
```

Note that for the Wyder/edgeR method, having only one replicate is also acceptable.

### Run:

Fill in `benchmark_imprinting/config/imprinting_config.txt` as per your requirements.

Run the following scripts: 
* `benchmark_imprinting/config/read_config_imprint.py benchmark_imprinting/config/imprinting_config.txt benchmark_imprinting/config/shell_env_imprint.txt`
* `benchmark_imprinting/wyder_imprinting.sh benchmark_imprinting/config/shell_env_imprint.txt`

### Output:

* `wyder_MEGs.txt` and `wyder_PEGs.txt` - imprinted genes called by Wyder/edgeR method
* `wyder_stats.txt` - edgeR stats for all genes

### References:
1.	Picard, C. L. & Gehring, M. Identification and Comparison of Imprinted Genes Across Plant Species. in Plant Epigenetics and Epigenomics 173–201 (Humana, New York, NY, 2020). doi:10.1007/978-1-0716-0179-2_13.

2.	Wyder, S., Raissig, M. T. & Grossniklaus, U. Consistent Reanalysis of Genome-wide Imprinting Studies in Plants Using Generalized Linear Models Increases Concordance across Datasets. Sci. Rep. 9, 1–13 (2019).

3.	Anderson, S. N., Zhou, P., Higgins, K., Brandvain, Y. & Springer, N. M. Widespread imprinting of transposable elements and young genes in the maize endosperm. Cold Spring Harbor Laboratory 2020.04.08.032573 (2020) doi:10.1101/2020.04.08.032573.


