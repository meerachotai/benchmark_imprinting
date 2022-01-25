#!/usr/bin/env Rscript

if (is.element('argparse', installed.packages()[,1])==FALSE) { stop("Error: R package argparse could not be loaded") }
suppressPackageStartupMessages(library(argparse))

parser = ArgumentParser()

parser$add_argument("outdir", nargs=1, help = "outdirectory")

args <- commandArgs(trailingOnly = TRUE)

opt <- parser$parse_args()

outdir = opt$outdir

true_megs = scan(paste0(outdir,"/true_MEGs.txt"), character())
true_pegs = scan(paste0(outdir,"/true_PEGs.txt"), character())
# all_genes = scan(paste0(outdir,"/all_genes.txt"), character())

if(file.exists(paste0(outdir,"/picard_MEGs.txt"))) {
	picard_megs = scan(paste0(outdir,"/picard_MEGs.txt"), character())
	picard_pegs = scan(paste0(outdir,"/picard_PEGs.txt"), character())
	
	true_picard_megs = intersect(true_megs, picard_megs)
	true_picard_pegs = intersect(true_pegs, picard_pegs)
	cat("-----------------------------------\n")
	cat("total Picard maternally-biased:", length(picard_megs), "\n")
	cat("true Picard maternally-biased: ",length(true_picard_megs), "\n") 
	cat("total Picard paternally-biased:", length(picard_pegs), "\n")
	cat("true Picard paternally-biased: ",length(true_picard_pegs), "\n")
	cat("-----------------------------------\n")
}

if(file.exists(paste0(outdir,"/wyder_MEGs.txt"))) {
	wyder_megs = scan(paste0(outdir,"/wyder_MEGs.txt"), character())
	wyder_pegs = scan(paste0(outdir,"/wyder_PEGs.txt"), character())
	
	true_wyder_megs = intersect(true_megs, wyder_megs)
	true_wyder_pegs = intersect(true_pegs, wyder_pegs)
	cat("-----------------------------------\n")
	cat("total Wyder maternally-biased:", length(wyder_megs), "\n")
	cat("true Wyder maternally-biased: ",length(true_wyder_megs), "\n") 
	cat("total Wyder paternally-biased:", length(wyder_pegs), "\n")
	cat("true Wyder paternally-biased: ",length(true_wyder_pegs), "\n")
	cat("-----------------------------------\n")
}

if(file.exists(paste0(outdir,"/anderson_MEGs.txt"))) {
	anderson_megs = scan(paste0(outdir,"/anderson_MEGs.txt"), character())
	anderson_pegs = scan(paste0(outdir,"/anderson_PEGs.txt"), character())
	
	true_anderson_megs = intersect(true_megs, anderson_megs)
	true_anderson_pegs = intersect(true_pegs, anderson_pegs)

	cat("-----------------------------------\n")
	cat("total Anderson maternally-biased:", length(anderson_megs), "\n")
	cat("true Anderson maternally-biased: ",length(true_anderson_megs), "\n") 
	cat("total Anderson paternally-biased:", length(anderson_pegs), "\n")
	cat("true Anderson paternally-biased: ",length(true_anderson_pegs), "\n")
	cat("-----------------------------------\n")
}

if(file.exists(paste0(outdir,"/roth_MEGs.txt"))) {
	roth_megs = scan(paste0(outdir,"/roth_MEGs.txt"), character())
	roth_pegs = scan(paste0(outdir,"/roth_PEGs.txt"), character())
	
	true_roth_megs = intersect(tolower(true_megs), megs)
	true_roth_pegs = intersect(tolower(true_pegs), pegs)

	cat("-----------------------------------\n")
	cat("total Roth maternally-biased:", length(roth_megs), "\n")
	cat("true Roth maternally-biased: ",length(true_roth_megs), "\n") 
	cat("total Roth paternally-biased:", length(roth_pegs), "\n")
	cat("true Roth paternally-biased: ",length(true_roth_pegs), "\n")
	cat("-----------------------------------\n")
}
