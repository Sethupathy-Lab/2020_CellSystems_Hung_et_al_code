#!/usr/bin/python3

import os
import sys
import re
import argparse
import textwrap
import errno
from collections import defaultdict

def check_file(file):
	"""
	check that file exists
	"""

	if not os.path.isfile(file):
		print("Error! File does not exist:", file)
		print("Exiting now...")
		return False

	return True


def gtf2bed(gtf, start, end, allgenes):
	"""
	convert gtf to bed file
	"""
	
	# check that files exist
	gtf_exist = check_file(gtf)
	if not gtf_exist:
		return

	# check that files have the right extension
	if not gtf.endswith(".gtf"):
		print("Error! GTF file extension doesn't appear right:", gtf)
		print("Exiting now...")
		return
	
	# read in gtf, loop through each line, for genes, print out new coordinates in bed6 format

	keep = ["protein_coding", "lincRNA", "antisense", "miRNA"]

	try:
		gtf_file = open(gtf)
	except:
		print("Error! Cannot open gtf file:", gtf)
		return

	# keep track of genes to check for duplicates
	genes = defaultdict(lambda: defaultdict(dict))
	genelog = []

	for line in gtf_file:
		if line.strip().startswith("#"): # skip comments
			continue

		chrom, source, feature, tss, stop, score, strand, frame, attribute = line.split('\t')[:9]
		
		if not feature == "gene":
			continue

		# make sure chromosome formatting is "chr1" and not "1"
		chrom = str(chrom)
		if chrom.isdigit() or chrom in ["X", "Y", "MT", "M"]:
			chrom = "chr" + chrom

		string1 = attribute.split('"; ')
		typepattern = re.compile("gene_.*type")
		typelist = list(filter(typepattern.search, string1))
		if len(typelist) != 1:
			print("Error! Maybe gtf attributes, especially gene type, are not formatted as expected!")
			print(typelist)
			return
		genetype = typelist[0].split('type "')[1].split('"')[0]

		ensemblpattern = re.compile("gene_id")
		ensembllist = list(filter(ensemblpattern.search, string1))
		if len(ensembllist) != 1:
			print("Error! Maybe gtf attributes, especially gene id, are not formatted as expected!")
			print(ensembllist)
			return
		ensembl = ensembllist[0].split('_id "')[1].split('"')[0]

		namepattern = re.compile("gene_name")
		namelist = list(filter(namepattern.search, string1))
		if len(namelist) != 1:
			print("Error! Maybe gtf attributes, especially gene name, are not formatted as expected!")
			print(namelist)
			return
		gene = namelist[0].split('name "')[1].split('"')[0]
		genestring = gene + "|" + ensembl + "|" + genetype

		# if not using all genes, filter genes
		if not allgenes: # if false
			if "pseudogene" not in str(genetype) and genetype not in keep:
				continue
		
		if strand == "+":
			tss = int(tss) - 1 + start # substract 1 b/c bed is 0-based and gtf is 1-based
			stop = int(stop) + end # don't subtract 1 b/c bed is half open and gtf is fully closed
		elif strand == "-":
			tss = int(tss) - 1 - end
			stop = int(stop) - start
		else:
			print("This gene is not on any strand:", gene)
			print("Exiting...")
			return

		if genestring in genelog:
			if chrom == 'chrX' or chrom == 'chrY': # x/y paralogs, let's keep both
				genestring = gene + "_" + chrom + "|" + genetype
				genes[chrom][genestring]['lower'] = tss
				genes[chrom][genestring]['upper'] = stop
				genes[chrom][genestring]['strand'] = strand
				genelog.append(genestring)
			elif genestring in genes[chrom]: # check if previously called gene is on the same chr
				if strand != genes[chrom][genestring]['strand']:
					sys.stderr.write("Error! Strand of this gene is not consistent between two entries: ")
					sys.stderr.write(genestring)
					sys.stderr.write("\n")
					continue
				genes[chrom][genestring]['lower'] = min(tss, genes[chrom][genestring]['lower'])
				genes[chrom][genestring]['upper'] = max(stop, genes[chrom][genestring]['upper'])
			else: # previously called gene is not on same chr! Error!
				sys.stderr.write("Error! Gene called multiple times and not on same chr: ")
				sys.stderr.write(genestring)
				sys.stderr.write(" ")
				sys.stderr.write(chrom)
				sys.stderr.write("\n")
				continue
		else:
			genes[chrom][genestring]['lower'] = tss
			genes[chrom][genestring]['upper'] = stop
			genes[chrom][genestring]['strand'] = strand
			genelog.append(genestring)

	gtf_file.close()		

	for key, value in genes.items():
		for key2 in value.keys():
			dist = genes[key][key2]['upper'] - genes[key][key2]['lower']
			lineout = key + "\t" + str(genes[key][key2]['lower']) + "\t" + str(genes[key][key2]['upper']) + "\t" + key2 + "\t" + str(dist) + "\t" + genes[key][key2]['strand'] + "\n"
			try:
				sys.stdout.write(lineout)
			except IOError as e: # handle python pipe error
				if e.errno == errno.EPIPE:
					pass

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Convert gtf file to bed file for gene coordinates. Output is in bed6 format. Output is printed to stdout so users can pipe to additional functions and/or redirect to a file.", formatter_class=argparse.RawDescriptionHelpFormatter, epilog=textwrap.dedent('''Common pipes and redirects: 
Look at first 20 genes					| head -n 20
Redirect to a file					> filename.bed'''))
	parser.add_argument("gtf", type=str, help="gencode gtf file")
	parser.add_argument("-s", "--start", type=int, help="distance in bp upstream (-) or downstream (+) of gene start to annotate (Default: 0 bp)", default=0)
	parser.add_argument("-e", "--end", type=int, help="distance in bp upstream (-) or downstream (+) of gene end to annotate (Default: 0 bp)", default=0)
	parser.add_argument("-a",  help="include all genes (Default: False. Will only include the following types of genes: protein coding, all types of pseudogenes, lincRNA, antisense, and miRNA)", action="store_true") 
	args = parser.parse_args()
	
	gtf2bed(args.gtf, args.start, args.end, args.a)
