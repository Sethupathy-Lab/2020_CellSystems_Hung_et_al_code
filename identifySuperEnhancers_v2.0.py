#!/usr/bin/python3

import os
import sys
import re
import argparse
import textwrap
import errno
import string
import random
import datetime
import subprocess
import shutil
from collections import defaultdict

def file_exist(files):
	"""
	check if each file in line-delimited list exists
	"""
	for bw in files:
		plus = bw.strip()
		minus = plus.replace('_plus.bw', '_minus.bw')
		
		if not os.path.isfile(plus):
			print("File doesn't exist: ", plus, flush = True)
			return False
		if not os.path.isfile(minus):
			print("File doesn't exist: ", minus, flush = True)
			return False
	
	return True


def read_gtf(gtf, up, down):
	"""
	read in gencode gtf file and promoter coordinates in dictionary
	"""
	transcripts = defaultdict(lambda: defaultdict(dict))

	print(datetime.datetime.now(), "Reading in gtf file...", flush = True)
	for line in gtf:
		if line.strip().startswith("#"): # skip comments
			continue

		chrom, source, feature, start, end, score, strand, frame, attribute = line.split('\t')
		
		if not feature == "transcript":
			continue

		# make sure chromosome formatting is "chr1" and not "1"
		chrom = str(chrom)
		if chrom.isdigit() or chrom in ["X", "Y", "MT", "M"]:
			chrom = "chr" + chrom

		string1 = attribute.split('"; ')
		transcriptpattern = re.compile("transcript_id")
		transcriptlist = list(filter(transcriptpattern.search, string1))
		if len(transcriptlist) != 1:
			print("Error! Maybe gtf attributes, especially transcript id, are not formatted as expected!")
			print(transcriptlist)
			return
		transcript = transcriptlist[0].split('id "')[1].split('"')[0]

		if strand == "+":
			transcripts[chrom][transcript]['lower'] = int(start) - up - 1 # subtract 1 b/c bed is 0-based and gtf is 1-based coordinates
			transcripts[chrom][transcript]['upper'] = int(start) + down - 1
		elif strand == "-":
			transcripts[chrom][transcript]['lower'] = int(end) - down # don't subtract one here b/c bed is half-open and gtf is fully closed
			transcripts[chrom][transcript]['upper'] = int(end) + up
		else:
			print("This transcript is not on any strand:", transcript, flush = True)
			print("Skipping", flush = True)
			continue
		
	gtf.close()
	return transcripts


def remove_proximal_TREs(transcripts, bed, tempbed): 
	"""
	removes proximal TREs from given bed file
	saves output as a tempoary bed file
	"""
	# read in bed file and only store distal (non-promoter) TREs
	
	# count total lines in bed file to keep track of progress
	counter = subprocess.run(["wc", "-l", bed.name], stdout = subprocess.PIPE).stdout.decode('utf-8').split(" ")[0]
	progress = [int(x*(int(counter))/20) for x in range(0, 21, 1)]
	progress_percent = list(range(0, 105, 5))
	
	# save distal TREs in a bed file
	unsorted = "unsorted_" + tempbed
	outfile = open(unsorted, 'w')

	linenum = 0
	for line in bed:
		# keep track of progress as script is running
		if linenum == progress[0]:
			percentage = str(progress_percent[0]) + "%"
			print(datetime.datetime.now(), "Removing proximal TREs...", percentage, "done", flush = True)
			progress.pop(0)
			progress_percent.pop(0)
			
		chromosome, first, last = line.split('\t')[:3]

		enhancer = True
		for key, value in transcripts[chromosome].items():
			
			if int(float(last)) <= int(transcripts[chromosome][key]['lower']):
				continue
			elif int(float(first)) >= int(transcripts[chromosome][key]['upper']):
				continue
			else:
				enhancer = False
				break
		
		# write only enhancers (distal TREs) to file
		if enhancer == True:
			outline = chromosome + "\t" + first + "\t" + last.rstrip() + "\n"
			outfile.write(outline)
		linenum += 1

	bed.close()
	outfile.close()

	command2 = "cat " + unsorted + " | sort -k1,1 -k2,2n > " + tempbed
	subprocess.run([command2], shell = True)
	os.remove(unsorted)


def stitch_enhancers(tempbed, stitch):
	"""
	stitch together enhancers within specified distance
	save stitched enhancers as temporary bed file
	determine enhancer lineage tree (i.e. which enhancer belongs to which stitched enhancer)
	output signal for each enhancer and indicate which stitched enhancer it belongs to
	"""

	# stitch together enhancers within given distance
	print(datetime.datetime.now(), "Stitching distal TREs together...", flush = True)
	tempbedstitched = tempbed[:57] + "_stitched.bed"
	stitchedfile = open(tempbedstitched, 'w')
	subprocess.run(["bedtools", "merge", "-i", tempbed, "-d", str(stitch)], stdout = stitchedfile)
	stitchedfile.close()

	# determine which enhancer belongs to which stitched enhancer
	counter = subprocess.run(["wc", "-l", tempbed], stdout = subprocess.PIPE).stdout.decode('utf-8').split(" ")[0]
	progress = [int(x*(int(counter))/20) for x in range(0, 21, 1)]
	progress_percent = list(range(0, 105, 5))

	# store enhancer counts in a dictionary
	counts = {}
	tempcounts = tempbed[:57] + "_counts.txt"
	countsfile = open(tempcounts, 'r')
	for line in countsfile:
		locus, signal = line.split("\t")
		counts[locus] = signal.strip()
	countsfile.close()

	# store each enhancer in a list
	single = open(tempbed, 'r')
	enhancer_list = []
	for line in single:
		enhancer_list.append(line)
	single.close()
	
	stitched = open(tempbedstitched, 'r')
	stitched_list = []
	for line in stitched:
		stitched_list.append(line)
	stitched.close()

	tempparentage = tempbed[:57] + "_parentage.txt"
	parentage = open(tempparentage, 'w')

	# loop through each enhancer and see it it belongs to each stitched enhancer
	parentage.write("TRE" + "\t" + "Signal" + "\t" + "Parent" + "\n")
	linenum = 0
	for line1 in enhancer_list:
		found = False
		if linenum == progress[0]:
			percentage = str(progress_percent[0]) + "%"
			print(datetime.datetime.now(), "Determining enhancer/stitched enhancer parentage...", percentage, "done", flush = True)
			progress.pop(0)
			progress_percent.pop(0)

		chr1, start1, stop1 = line1.split('\t')[:3]
		outline1 = chr1 + ":" + start1 + "-" + stop1.strip()
		for line2 in stitched_list:
			chr2, start2, stop2 = line2.split('\t')[:3]
			
			if chr1 != chr2:
				continue
			elif int(float(start1)) >= int(start2) and int(float(stop1)) <= int(stop2):
				if found == True:
					print("Error! Single enhancer belongs to multiple stitched enhancers:", outline1, flush = True)
					print("Exiting...", flush = True)
					return
				else:
					found == True
					parent = "\t" + chr2 + ":" + start2 + "-" + stop2.strip() + "\n"
					fullout = outline1 + "\t" + counts[outline1] + parent
					parentage.write(fullout)
			else:
				continue

		linenum += 1
	parentage.close()




def super_enhancers(bed, gtf, bigwig, up, down, stitch, alltres, prefix):
	"""
	identify super enhancers from ChRO-seq data
	"""

	starttime = datetime.datetime.now()
	
	print(flush = True)	
	print("TRE bed file:", bed.name, flush = True)
	print("GTF file:", gtf.name, flush = True)
	print("File with list of bigwigs:", bigwig.name, flush = True)
	print("Keep all TREs:", alltres, flush = True)
	print("Remove proximal TREs:", not alltres, flush = True)
	print("Proximal TREs extend upstream of TSS:", up, "bp", flush = True)
	print("Proximal TREs extend downstream of TSS:", down, "bp", flush = True)
	print("TRE stitching distance:", stitch, "bp", flush = True)
	print("Output file prefix given:", prefix, flush = True)
	print(flush = True)

	exists = file_exist(bigwig)
	if not exists:
		print("Exiting...", flush = True)
		return
	
	tempbed = ''.join(random.choices(string.ascii_letters + string.digits, k = 48)) + '_enhancer.bed'
	if prefix == 'FALSE':
		prefix = tempbed[:48]

	file1 = prefix + "_stitched_enhancers_info.txt"
	file2 = prefix + "_SuperEnhancersRanked.png"
	if os.path.isfile(file1) or os.path.isfile(file2):
		print("File outputs already exist! Please don't make me overwrite them! Check for:", file1, "and/or", file2, flush = True)
		print("Exiting...", flush = True)
		return

	print("Input plus strand bigwig files:", flush = True)
	bigwigfiles = open(bigwig.name, 'r')
	linetracker = 1
	for bw in bigwigfiles:
		number = str(linetracker) + ":"
		print(number, bw.strip(), flush = True)
		linetracker += 1
	print(flush = True)
	bigwigfiles.close()

	if not alltres:
		transcripts = read_gtf(gtf, up, down)
		remove_proximal_TREs(transcripts, bed, tempbed)
	else:
		print(datetime.datetime.now(), "Not removing proximal TREs...", flush = True)
		shutil.copyfile(bed.name, tempbed)
	
	print(datetime.datetime.now(), "Calling R script to extract read counts at distal TREs...", flush = True)
	command = "R --vanilla --slave --args " + bigwig.name + " " + tempbed + " < /home/pr46_0001/ChROseq/ChROseq_pipeline/Tools/bin/getCounts.R"
	subprocess.run([command], shell = True)

	stitch_enhancers(tempbed, stitch)

	print(datetime.datetime.now(), "Calling R script to call super enhancers...", flush = True)
	command2 = "R --vanilla --slave --args " + tempbed + " " + prefix +  " < /home/pr46_0001/ChROseq/ChROseq_pipeline/Tools/bin/SuperEnhancers.R"
	subprocess.run([command2], shell = True)

	print(datetime.datetime.now(), "Removing temporary files...", flush = True)
	os.remove(tempbed)
	os.remove(tempbed[:57] + "_stitched.bed")
	os.remove(tempbed[:57] + "_counts.txt")
	os.remove(tempbed[:57] + "_parentage.txt")

	endtime = datetime.datetime.now()
	elapsed = endtime - starttime
	elapsed = divmod(elapsed.total_seconds(), 60)
	print(datetime.datetime.now(), "Done! Runtime:", int(elapsed[0]), "min,", elapsed[1], "seconds", flush = True)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''Identify Super/Stretch enhancers from ChRO-seq data. This script identifies Super/Stretch enhancers from a single condition by doing the following things:
1. Removing all proximal TREs (optional)
2. Stitch together TREs within defined distance
3. Get read counts within each TRE from bigwigs
4. Normalizing read counts for differences in sequencing depth
5. Sum TRE read counts within each stitched enhancer
6. Rank enhancers and identify super enhancers''')) 
	parser.add_argument("bed", type=argparse.FileType('r'), help="bed file of TREs. Accepts a file or stdin", nargs="?", default=sys.stdin)
	parser.add_argument("gtf", type=argparse.FileType('r'), help="gencode gtf file")
	parser.add_argument("bigwig", type=argparse.FileType('r'), help="line-delimited file with path to bigwig file for each sample of interest. Must include absolute (full) path to each bigwig file. Include only path to plus strand bigwig. Minus strand bigwig must have same prefix")
	parser.add_argument("-u", "--up", type=int, help="distance in bp upstream of TSS to consider promoter (Default: 1000 bp)", default=1000)
	parser.add_argument("-d", "--down", type=int, help="distance in bp downstream of TSS to consider promoter (Default: 0 bp)", default=0)
	parser.add_argument("-s", "--stitch", type=int, help="distance in bp to use as TRE stiching distance. TREs within this distance will be stiched together. For Super enhancers, 12500 bp is typical. For Stretch enhancers, 3000 bp is typical, although defining Stretch enhancers does not use enhancer stitching, but rather a length threshold of >= 3 kb of enhancer-defined states from ChromHMM (Default: 12500 bp)", default=12500)
	parser.add_argument("-p", "--prefix", type=str, default='FALSE', help="prefix for output files (Default: random string)")
	parser.add_argument("-a", "--all", action="store_true", help="include all TREs in analysis. Ignores promoter definition defined by -u/-d. Not recommended for run-on sequencing data due to presence of promoter proximal pausing not related to enhancers")
	args = parser.parse_args()
	
	super_enhancers(args.bed, args.gtf, args.bigwig, args.up, args.down, args.stitch, args.all, args.prefix)
