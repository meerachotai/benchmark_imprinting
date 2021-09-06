#!/usr/bin/env python

import random
import scipy
from scipy.stats import powerlaw
import argparse
import numpy

parser = argparse.ArgumentParser()
parser.add_argument('-i', default = 3, type=float,help = 'indel score')
parser.add_argument('-s', default = 1, type=float, help = 'snp score')
parser.add_argument('-m', default = 1, type=float, help = 'match score')
parser.add_argument('-e', default = 2, type=float,help = 'extend score')
parser.add_argument('-r', default = 1, type=int, help = 'snp:indel ratio')
parser.add_argument('-S', default = 97, type=float, help = 'score')
parser.add_argument('-d', default = 3, type=int, help = "standard dev. for similarity score dist.")
parser.add_argument('-w', default = 5, type=float, help = 'indel score window')
parser.add_argument('-n', type=int, help = 'number of chromosomes')
parser.add_argument('-a', default = 0.2, type=float, help = 'indel length a value')
parser.add_argument('-R', default = None, type=int, help = "seed")
parser.add_argument('-A', default = "", help = 'strain A FASTA file name')
parser.add_argument('-B', default = "", help = 'strain B FASTA file name')
parser.add_argument('-o', default = "", help = 'outprefix')

args = parser.parse_args()

match = args.m
snp = args.s
indel = args.i
extend = args.e
si_ratio = args.r
mean_score = args.S
window=args.w
seqA = args.A
standard_dev = args.d
n = args.n
alpha = args.a
seed = args.R

seqB = args.B
ref2seq = args.o + "_ref2sim.txt"
similarity = args.o + "_similarity.txt"


print("Number of chromosomes:", n)
print("Achieving mean similarity:", mean_score, "%")
print("with normal dist. standard dev. of:", standard_dev)
print("Score settings:")
print("match:", match)
print("SNP:", snp)
print("indel:", indel, "extend/gap:", extend)
print("SNP:indel ratio of:", si_ratio)
print("With an indel score window of:", window)
print("Indel length under the power-law dist. = a*x^(a-1) where a =", alpha)
print("Reading file:", seqA)
print("Writing to files:", seqB, ref2seq, similarity)
if(seed != None):
	print("For reproducibility, using seed:", seed)
else:
	print("using current system time as default seed (cannot reproduce)")

###################
#### functions ####
###################


def add_snps(snp_count, sequence, length):
    
	random.seed(a = seed) # set seed
	bases = ["A", "T", "G", "C"]
	snps = []
	pos_chosen = []

	for i in range(snp_count):
		ref_pos = random.choice([i for i in range(1,length) if i not in pos_chosen])
		pos_chosen.append(ref_pos)
		ref_allele = sequence[ref_pos-1]
		pos = random.randint(0,2)
		sim_allele = [base for base in bases if base != ref_allele][pos]

		added_snp = [ref_pos, ref_pos, ref_allele, ref_pos, ref_pos, sim_allele, "SNP",0]
		snps.append(added_snp)
        
	return snps

def add_indels(need_indel_score, sequence, length, window):
    
    random.seed(a = seed)  
#     scipy_randomGen.random_state=seed
    
    bases = ["A", "T", "G", "C"]
    indels = []
    indel_score = 0
    
    while(indel_score < need_indel_score):
        
        indel_length = scipy.stats.powerlaw.rvs(a = alpha, loc = 1, scale = 9, size = 1, random_state = numpy_randomGen) # range - 1 - 10
        indel_length = round(indel_length[0])
        indel_length = int(indel_length)
        new_indel_score = indel_score + ((indel_length - 1) * extend) + indel 
        
        while(new_indel_score > need_indel_score + window): # so the last one isn't WAY past the needed score
            indel_length = scipy.stats.powerlaw.rvs(a = alpha, loc = 1, scale = 9, size = 1, random_state = numpy_randomGen)
            indel_length = round(indel_length[0])
            indel_length = int(indel_length)
            new_indel_score = indel_score + ((indel_length - 1) * extend) + indel 
            
        indel_score = new_indel_score
            
        
        insertion = random.randint(0,1) # 0 = deletion, 1 = insertion
        
        if(insertion):
            ref_start = random.randint(1,length)
            ref_end = ref_start
            ref_allele = sequence[ref_start-1]
            
            sim_allele = ref_allele
            sim_start = ref_start
            sim_end = sim_start + indel_length
            for i in range(indel_length):
                pos = random.randint(0,3)
                sim_allele += bases[pos]
    
            change = indel_length # must track this while changing new sequence !!!!!!
        
        else:
            ref_start = random.randint(1,length-indel_length)
            ref_end = ref_start+indel_length
            ref_allele = sequence[ref_start-1:ref_end] # 0-indexed, so subtract 1, NOT INCLUDING last index
            
            sim_allele = ref_allele[0]
            sim_start = ref_start
            sim_end = ref_start
            
            change = -indel_length
        
        
        added_indel = [ref_start, ref_end, ref_allele, sim_start, sim_end, sim_allele, "INDEL", change]
        indels.append(added_indel) 
        
    return indels, indel_score

# new_seq = new_seq[0:ref_start-1+pos] + sim_allele + new_seq[ref_start+pos:]
def ref_to_sim(snps, indels, chrom):   
    pos = 0
    new_seq = chrom

    variants = snps + indels
    variants.sort()

    for variant in variants:

        ref_start = variant[0]
        change = variant[-1]
        if(variant[-2] == "SNP"):
            sim_allele = variant[5]
            new_seq = new_seq[0:ref_start-1+pos] + sim_allele + new_seq[ref_start+pos:]
        else:
            if(change < 0): # deletion
            	sim_allele = variant[5]
            	ref_end = variant[1]
            	new_seq = new_seq[0:ref_start+pos] + new_seq[ref_end+pos:]
            	variant[5] = new_seq[ref_start-1+pos] + sim_allele[1:]
            else: # insertion
            	sim_allele = variant[5]
            	new_seq = new_seq[0:ref_start+pos] + sim_allele[1:] + new_seq[ref_start+pos:]
            	variant[5] = new_seq[ref_start-1+pos] + sim_allele[1:] # in case there's a snp in that pos, sync
        # change sim_start and sim_end
        variant[3] += pos
        variant[4] += pos

        pos += change

    return new_seq, variants

def edit_genome(header, chrom, write_seq, write_txt, write_sim, score_percent, window):
        
#         print(header, chrom)
        length = len(chrom)
        perfect_score = length * match
        score = round((score_percent/100) * perfect_score)
        print("--------------------------------")
        print("Chromosome", counter)
        print("length of chromosome:", length)
        print("perfect score:", perfect_score)
        print("required score:", score_percent, "% or", score, "/", perfect_score)
        
        snp_count = round((perfect_score - score) / (match + snp + (snp/si_ratio)))
        snp_count = int(snp_count)
        snp_score = snp_count * snp
        need_indel_score = round((snp_count * snp) / si_ratio)
        print("number of SNPs needed:", snp_count)
        print("SNP score needed:", snp_score)
        print("indel score needed:", need_indel_score)
        
        indels, indel_score = add_indels(need_indel_score,chrom, length, window)
        print("actual indel score:", indel_score)
        actual_score = ((length - snp_count) * match) - snp_score - indel_score
        actual_percent = (actual_score * 100) / perfect_score
        print("actual score:", actual_percent, "% or", actual_score, "/", perfect_score)
        snps = add_snps(snp_count, chrom, length)
        new_seq, variants = ref_to_sim(snps, indels, chrom)
        
        header = header.split(" ")[0].strip()
        write_seq.write(header)
        write_seq.write("\n")
        write_seq.write(new_seq)
        write_seq.write("\n")
        
        header = header[1:].strip() + "\t"
        for variant in variants:
            variant = variant[:-1]
            str_variant = [str(x) for x in variant]
            write_txt.write("\n")
            write_txt.write(header)
            write_txt.write('\t'.join(str_variant[0:3]))
            write_txt.write("\t")
            write_txt.write(header)
            write_txt.write('\t'.join(str_variant[3:]))
        
        write_sim.write(header)
        write_sim.write(str(actual_percent))
        write_sim.write("\n")

def conserved_genome(header, chrom, write_seq, write_sim):
    print("--------------------------------")
    print("Chromosome", counter, ": conserved!")

    header = header.split(" ")[0].strip()
    write_seq.write(header)
    write_seq.write("\n")
    write_seq.write(chrom)
    write_seq.write("\n")

    header = header[1:].strip() + "\t"
    write_sim.write(header)
    write_sim.write("100")
    write_sim.write("\n")

read_seq = open(seqA, "r")
write_seq = open(seqB, "w")
write_txt = open(ref2seq, "w")
write_sim = open(similarity, "w")

txt_header = ["ref_chr", "ref_start", "ref_end", "ref_allele", "sim_chr", "sim_start", "sim_end", "sim_allele", "variant"]
write_txt.write("\t".join(txt_header))        
write_sim.write("chrom\tsimilarity%\n")

chrom = ""
header = read_seq.readline()


numpy_randomGen = numpy.random.RandomState(seed=seed) # set seed
scores_norm = scipy.stats.norm.rvs(loc = mean_score, scale = standard_dev, size = n, random_state = numpy_randomGen)

counter = 1
for line in read_seq:
    if(">" not in line):
        line = line.strip()
        chrom += line
    else:
        req_score = scores_norm[counter-1]
        if(req_score < 100):
            edit_genome(header, chrom, write_seq, write_txt, write_sim,req_score, window)
        else: # no change, conserved
            conserved_genome(header, chrom, write_seq, write_sim)
        # refresh for next chrom
       
        counter += 1
        header = line
        chrom = ""
        
# for the last chrom
req_score = scores_norm[counter-1]
if(req_score < 100):
    edit_genome(header, chrom, write_seq, write_txt, write_sim, req_score, window)
else: # no change, conserved
    conserved_genome(header, chrom, write_seq, write_sim)
         
            
read_seq.close()
write_seq.close()
write_txt.close()
write_sim.close()

print("Done!")
