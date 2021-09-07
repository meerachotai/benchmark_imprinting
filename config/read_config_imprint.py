#!/usr/bin/env python3

import random
import sys

input_file = str(sys.argv[1])
output_env = str(sys.argv[2])

# input_file = "imprinting_config.txt"
inputs = []
complete = True

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
				next_in = next_in.lower()
			inputs.append(next_in)

file = open(output_env, "w")

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
file.write("mat_cutoff_picard="+inputs[10]+"\n")
file.write("pat_cutoff_picard="+inputs[11]+"\n")
file.write("paired="+inputs[12]+"\n")
file.write("rep="+inputs[13]+"\n")
file.write("AxB_rep="+inputs[14]+"\n")
file.write("BxA_rep="+inputs[15]+"\n")
file.write("majority="+inputs[16]+"\n")

file.write("counts_dir="+inputs[17]+"\n")
file.write("pval_wyder="+inputs[18]+"\n")
file.write("logfc_wyder="+inputs[19]+"\n")

file.write("refA="+inputs[20]+"\n")
file.write("annotA="+inputs[21]+"\n")
file.write("refB="+inputs[22]+"\n")
file.write("annotB="+inputs[23]+"\n")
file.write("gene_key="+inputs[24]+"\n")
file.write("htseq_i="+inputs[25]+"\n")
file.write("pval="+inputs[26]+"\n")
file.write("mat_cutoff_anderson="+inputs[27]+"\n")
file.write("pat_cutoff_anderson="+inputs[28]+"\n")
file.write("logfc="+inputs[29]+"\n")
file.write("rep="+inputs[30]+"\n")
file.write("rename="+inputs[31]+"\n")
file.write("a_annot="+inputs[32]+"\n")
file.write("b_annot="+inputs[33]+"\n")

file.close()
