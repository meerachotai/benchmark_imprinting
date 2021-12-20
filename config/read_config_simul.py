#!/usr/bin/env python3

import random
import sys

input_file = str(sys.argv[1])
output_env = str(sys.argv[2])

inputs = []
# input_file = "simulation_config.txt"

with open(input_file, 'r') as f:
     
	for line in f:
		if line.startswith("#"):
			continue
		elif line.startswith("\n") or line.isspace():
			continue
		else:
			split_line = line.split('=')
			next_in = split_line[1].strip()
		if(next_in.isspace() or next_in == ''):
			if(split_line[0].strip() == "SEED"): # for simulating
				seed = random.randint(1, 1234567)
				inputs.append(seed)
				print("Using randomly-generated seed", seed)
			else:
				print(split_line[0], "is empty, fill in config file and try again")
				quit()
		if(split_line[0].strip() == "DISPERSION"):
			next_in = next_in.lower()
		inputs.append(next_in)

file = open(output_env, "w")

file.write("outdir="+inputs[0] + "\n")
file.write("scripts_dir="+inputs[1]+"\n")
file.write("total_n="+inputs[2]+"\n")
file.write("ref="+inputs[3]+"\n")
file.write("annot="+inputs[4]+"\n")
file.write("strainA="+inputs[5]+"\n")
file.write("strainB="+inputs[6]+"\n")
file.write("refA="+inputs[7]+"\n")
file.write("refB="+inputs[8]+"\n")
file.write("annotA="+inputs[9]+"\n")
file.write("annotB="+inputs[10]+"\n")
file.write("seed="+inputs[11]+"\n")
file.write("sim_score="+inputs[12]+"\n")
file.write("window="+inputs[13]+"\n")
file.write("snp_score="+inputs[14]+"\n")
file.write("indel_score="+inputs[15]+"\n")
file.write("extend_score="+inputs[16]+"\n")
file.write("match_score="+inputs[17]+"\n")
file.write("ratio="+inputs[18]+"\n")
file.write("unbiased="+inputs[19]+"\n")
file.write("meg="+inputs[20]+"\n")
file.write("peg="+inputs[21]+"\n")
file.write("megs_bias="+inputs[22]+"\n")
file.write("pegs_bias="+inputs[23]+"\n")
file.write("disp="+inputs[24]+"\n")
file.write("rep="+inputs[25]+"\n")
file.write("read_length="+inputs[26]+"\n")


file.close()
