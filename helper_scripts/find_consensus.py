#!/usr/bin/env python3

import argparse, math

parser = argparse.ArgumentParser()
parser.add_argument('-i', default = "", help = 'input all MEGs/PEGs file name')
parser.add_argument('-t', default = 0, type=int, help = 'total reciprocal cross pairs')
parser.add_argument('-m', default = 0, type=int, help = 'majority reciprocal cross pairs needed for consensus')
parser.add_argument('-o', default = "", help = 'output file with consensus MEGs/PEGs list')

args = parser.parse_args()

input_file = args.i
total = args.t
majority = args.m
output_file = args.o

if(total != 0):
	majority = math.ceil(total/2)
	print("Using majority voting of >= ", majority, " for consensus calls")
	
imprinted = {}
consensus = []

print("reading input MEGs/PEGs list...")
with open(input_file, 'r') as rFile:
	for line in rFile:
		line = line.split("\t")[0]
		if line not in imprinted:
			imprinted[line] = 1
		else:
			imprinted[line] += 1
		

for gene in imprinted:
	if imprinted[gene] >= majority:
		consensus.append(gene)

print("writing output consensus list...")
with open(output_file, 'w') as wFile:
	wFile.write('\n'.join(consensus))
