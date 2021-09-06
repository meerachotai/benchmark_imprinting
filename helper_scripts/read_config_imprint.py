#!/usr/bin/env python3

# if(split_line[0].strip() == "SNP_FILE" or split_line[0].strip() == "FASTQ_DIR" or split_line[0].strip() == "AxB_REP" or split_line[0].strip() == "BxA_REP" or split_line[0].strip() == "REP"): # for PICARD
# 	inputs.append("")

import random
inputs = []
complete = True
input_file = "imprinting_config.txt"

with open(input_file, 'r') as f:
	for line in f:
		if line.startswith("#"):
			continue
		elif line.startswith("\n") or line.isspace():
			continue
		else:
			split_line = line.split('=')
# 			print(split_line)
			next_in = split_line[-1].strip()
			if(next_in.isspace()):
				next_in = ""
			if(next_in.lower() == "true" or next_in.lower() == "false"):
				next_in.lower()
			inputs.append(next_in)

file = open("shell_env_imprint.txt", "w")

file.write("outdir="+inputs[0] + "\n")
file.write("scripts_dir="+inputs[1]+"\n")

file.write("strainA="+inputs[2]+"\n")
file.write("strainB="+inputs[3]+"\n")
file.write("fastq_dir="+inputs[4]+"\n")
file.write("paired_end="+inputs[5]+"\n")

file.write("picard="+inputs[6]+"\n")
file.write("genome="+inputs[7]+"\n")
file.write("annot="+inputs[8]+"\n")
file.write("snps="+inputs[9]+"\n")
file.write("mat_cutoff="+inputs[10]+"\n")
file.write("pat_cutoff="+inputs[11]+"\n")
file.write("paired="+inputs[12]+"\n")
file.write("rep="+inputs[13]+"\n")
file.write("AxB_rep="+inputs[14]+"\n")
file.write("BxA_rep="+inputs[15]+"\n")
file.write("majority="+inputs[16]+"\n")

file.write("refA="+inputs[17]+"\n")
file.write("annotA="+inputs[18]+"\n")
file.write("refB="+inputs[19]+"\n")
file.write("annotB="+inputs[20]+"\n")
file.write("gene_key="+inputs[21]+"\n")
file.write("htseq_i="+inputs[22]+"\n")
file.write("pval="+inputs[23]+"\n")
file.write("mat_cutoff="+inputs[24]+"\n")
file.write("pat_cutoff="+inputs[25]+"\n")
file.write("logfc="+inputs[26]+"\n")
file.write("rep="+inputs[27]+"\n")
file.write("rename="+inputs[28]+"\n")
file.write("a_annot="+inputs[29]+"\n")
file.write("b_annot="+inputs[30]+"\n")

file.close()