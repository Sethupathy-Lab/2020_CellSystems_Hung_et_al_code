#!/usr/bin/python3

import os
import sys
import argparse
import textwrap
import errno
import string
import random
import subprocess
import re

def check_file(file):
	"""
	check that file exists
	"""

	if not os.path.isfile(file):
		print("Error! File does not exist:", file)
		print("Exiting now...")
		return False

	return True


def closest_gene(bed, gtf, exp, useGene):
	"""
	find closest gene to each TRE
	"""

	# check that files exist
	bed_exist = check_file(bed)
	gtf_exist = check_file(gtf)
	if not (bed_exist and gtf_exist):
		return

	if exp != None:
		exp_exist = check_file(exp)
		if not (exp_exist):
			return

	# check that files have the right extension
	if not bed.endswith(".bed"):
		print("Error! Bed file extension doesn't appear right:", bed)
		print("Exiting now...")
		return
	if not gtf.endswith(".gtf"):
		print("Error! GTF file extension doesn't appear right:", gtf)
		print("Exiting now...")
		return

	# read in expressed genes
	if exp != None:
		try:
			exp_genes = open(exp)
		except:
			print("Error! Cannot open list of expressed genes:", exp)
			return

		expressed = list()
		for line in exp_genes:
			exp_gene_name = line.strip().upper()
			expressed.append(exp_gene_name)


	# read in gtf and store gene coordinates and gene info as bed file

	try:
		gtf_file = open(gtf)
	except:
		print("Error! Cannot open gtf file:", gtf)
		return

	# make temporary bed file
	makebed = True
	while makebed == True:
		tempbed = ''.join(random.choices(string.ascii_letters + string.digits, k = 48)) + '.bed'
		if not os.path.isfile(tempbed):
			makebed = False
	tempbedsorted = tempbed[:-4] + ".sorted.bed"

	tempout = open(tempbed, 'w')

	for line in gtf_file:
		if line.strip().startswith("#"): # skip comments
			continue

		chrom, source, feature, start, end, score, strand, frame, attribute = line.split('\t')
		
		if not feature == "gene":
			continue

		# make sure chromosome formatting is "chr1" rather than "1" to match with bed formatting
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

		namepattern = re.compile("gene_name")
		namelist = list(filter(namepattern.search, string1))
		if len(namelist) != 1:
			print("Error! Maybe gtf attributes, especially gene name, are not formatted as expected!")
			print(namelist)
			return
		gene = namelist[0].split('name "')[1].split('"')[0]

		if exp != None: # if using only expressed gene list
			if not gene.strip().upper() in expressed: # if gene is not in expressed list
				continue
		if useGene: # use entire gene for closest calculation
			start1 = str(int(start) - 1) # substract 1 b/c gtf is 1-based and bed is 0-based. don't need to do this for end b/c gtf is fully closed (inclusive) and bed is half open (end not included)
			end1 = end
		else: # use TSS for closest calculation
			if strand == "+":
				start1 = str(int(start) - 1)
				end1 = start
			else:
				start1 = str(int(end) - 1)
				end1 = end

		outline = '\t'.join([chrom, start1, end, gene, genetype, strand])
		tempout.write(outline)
		tempout.write("\n")

	tempout.close()
	gtf_file.close()

	# sort new temporary bed file
	sortedbed = subprocess.run(["sort", "-k1,1", "-k2,2n", tempbed], stdout = subprocess.PIPE)
	
	tempsortout = open(tempbedsorted, 'w')
	tempsortout.write(sortedbed.stdout.decode('utf-8'))
	tempsortout.close()

	# make sure TRE bed file is sorted
	tresorted = tempbed[:48] + "_tre.bed"
	command2 = "cat " + bed + " | sort -k1,1 -k2,2n > " + tresorted
	subprocess.run([command2], shell = True)

	closestgene = subprocess.run(["bedtools", "closest", "-a", tresorted, "-b", tempbedsorted, "-D", "b"], stdout = subprocess.PIPE)
	
	os.remove(tempbed)
	os.remove(tempbedsorted)

	# get number of rows in bed file
	linenumber = 1
	try:
 		bed_file = open(tresorted, 'r')
	except:
 		print("Error! Cannot open gtf file:", gtf)
 		return

	for line in bed_file:
		if linenumber == 10: # doesn't really matter what row, just go down a few in case of comments or something weird
			numcolumns = len(line.split('\t'))
		if linenumber > 10:
			break
		linenumber += 1
	

	previousTRE = ""
	previousLine = ""
	tracker = False
	#bedtoolsout = closestgene.stdout.decode('utf-8').splitlines()
	for index, line in enumerate(closestgene.stdout.decode('utf-8').splitlines(), start = 1):
		line_elements = line.split('\t')
		closestcolnum = len(line_elements)
		tre = ''
		for x in range(numcolumns):
			tre = '.'.join([tre, line_elements[x]])
		tre = tre.strip('.')
		# get direction of TRE relative to gene
		if int(line_elements[closestcolnum - 1]) == 0:
			direction = "Overlapping"
		elif int(line_elements[closestcolnum - 1]) < 0:
			direction = "Upstream"
		else: # int(line_elements[closestcolnum - 1]) > 0:
			direction = "Downstream"
		
		# determine if TRE is equally close to or overlapping multiple genes
		if tracker == True:
			previousLine = '\t'.join([previousLine, "Tie"])
			if tre == previousTRE:
				tracker = True
			else:
				tracker = False
		else:
			if tre == previousTRE:
				previousLine = '\t'.join([previousLine, "Tie"])
				tracker = True
			else:
				previousLine = '\t'.join([previousLine, "-"])
				tracker = False
		previousTRE = tre

		if index > 1:
			try:
				sys.stdout.write(previousLine + '\n')
			except IOError as e: # handle python pipe error
				if e.errno == errno.EPIPE:
					pass	

		previousLine = '\t'.join([line, direction])
	
		if index == len(closestgene.stdout.decode('utf-8').splitlines()):
			if tracker == True:
				previousLine = '\t'.join([previousLine, "Tie"])
			else:
				previousLine = '\t'.join([previousLine, "-"])
			try:
				sys.stdout.write(previousLine + '\n')
			except IOError as e: # handle python pipe error
				if e.errno == errno.EPIPE:
					pass
			continue

	bed_file.close()
	os.remove(tresorted)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Find closest gene to each TRE based on TSS. For ties (multiple genes that are equally close to or overlap a TRE, all tied genes are listed. Output includes the TRE bed file followed by gene information including gene type, strand, distance from gene, and whether or not there are other tied genes for that TRE. Output is printed to stdout so users can pipe to additional functions and/or redirect to a file.", formatter_class=argparse.RawDescriptionHelpFormatter, epilog=textwrap.dedent('''Common pipes and redirects: 
Sort by genes with close TREs from high to low		| awk "{print $9}' | sort | uniq -c | sort -k1,1nr
Extract all lincRNA (or other class) genes		| awk '$10 == "lincRNA" {print $0}'
Look at first 20 genes					| head -n 20
Redirect to a file					> filename.txt'''))
	parser.add_argument("bed", type=str, help="bed file of TREs")
	parser.add_argument("gtf", type=str, help="gencode gtf file")
	parser.add_argument("-e", "--express", type=str, help="consider only expressed genes. Include a line-delimited file of expressed genes to use")
	parser.add_argument("-g", "--gene", action="store_true", help="find closest gene based on entire gene--either TSS or end. (Default: use TSS only)")
	args = parser.parse_args()
	
	closest_gene(args.bed, args.gtf, args.express, args.gene)
