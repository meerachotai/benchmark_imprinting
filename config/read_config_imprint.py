#!/usr/bin/env python3

import random
import sys
import pandas as pd

input_file = str(sys.argv[1])
output_env = str(sys.argv[2])

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
file.write("stranded="+inputs[6]+"\n")

if(inputs[4] == ""):
	if(inputs[7] == ""):
		print("Please provide either a FASTQ directory prefix under FASTQ_DIR or provide a tab-delimited config file under FASTQ_CONFIG_FILE")
		quit()
	else:
		config = pd.read_table(inputs[7],delimiter = "\t", header=None)
		config = config.replace({config.columns[2]:'AxB'},1)
		config = config.replace({config.columns[2]:'BxA'},2)
		config = config.sort_values(by=[config.columns[3]], axis = 0)
		config = config.sort_values(by=[config.columns[1]], axis = 0)
		config = config.replace({config.columns[2]:'F'},1)
		config = config.replace({config.columns[2]:'R'},2)
		config = config.reset_index(drop=True)	

		strand = config[config.columns[2]].tolist()
		if(2 in strand): # reverse file is given
			file.write("paired_end=true\n")
			for i,s in enumerate(strand):
				if(i % 2):
					if(s == 1):
						b, c = config.iloc[i-1].copy(), config.iloc[i].copy()
						config.iloc[i],config.iloc[i-1] = b,c
		filenames = config[[config.columns[0]]]
		filenames.to_csv("rearranged_config.txt",sep='\t', index=False, header=False)
		file.write("config=rearranged_config.txt\n")
        
file.write("picard="+inputs[8]+"\n")
file.write("genome="+inputs[9]+"\n")
file.write("annot="+inputs[10]+"\n")
file.write("snps="+inputs[11]+"\n")
file.write("mat_cutoff_picard="+inputs[12]+"\n")
file.write("pat_cutoff_picard="+inputs[13]+"\n")
file.write("pval_picard="+inputs[14]+"\n")
file.write("paired="+inputs[15]+"\n")
file.write("rep="+inputs[16]+"\n")
file.write("AxB_rep="+inputs[17]+"\n")
file.write("BxA_rep="+inputs[18]+"\n")
file.write("majority="+inputs[19]+"\n")

file.write("counts_dir_wyder="+inputs[20]+"\n")
file.write("pval_wyder="+inputs[21]+"\n")
file.write("logfc_wyder="+inputs[22]+"\n")

file.write("refA="+inputs[23]+"\n")
file.write("annotA="+inputs[24]+"\n")
file.write("refB="+inputs[25]+"\n")
file.write("annotB="+inputs[26]+"\n")
file.write("gene_key="+inputs[27]+"\n")
file.write("htseq_i="+inputs[28]+"\n")
file.write("pval_anderson="+inputs[29]+"\n")
file.write("mat_cutoff_anderson="+inputs[30]+"\n")
file.write("pat_cutoff_anderson="+inputs[31]+"\n")
file.write("logfc_anderson="+inputs[32]+"\n")
file.write("rep="+inputs[33]+"\n")
file.write("rename="+inputs[34]+"\n")
file.write("a_annot="+inputs[35]+"\n")
file.write("b_annot="+inputs[36]+"\n")

file.write("counts_dir_roth="+inputs[37]+"\n")
file.write("pval_roth="+inputs[38]+"\n")
file.write("logfc_roth="+inputs[39]+"\n")
file.write("cutoff="+inputs[40]+"\n")

file.close()
