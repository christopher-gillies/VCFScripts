#!/usr/bin/python

import argparse
import re
import os
import string
from pysam import VariantFile
import itertools

def four_gamete_rule(counts):
	n = sum(counts)
	num_haps = float(2 * n)
	n00,n01,n02,n10,n11,n12,n20,n21,n22 = counts
	freq_hap_00 =  float(2 * n00 + n10 + n01) / num_haps
	freq_hap_01 =  float(2 * n02 + n01 + n12) / num_haps
	freq_hap_10 =  float(2 * n20 + n21 + n10) / num_haps
	freq_hap_11 =  float(2 * n22 + n21 + n12) / num_haps
	hap_freqs = (freq_hap_00,freq_hap_01,freq_hap_10,freq_hap_11)
	return min(hap_freqs)

def four_gamete_rule_test():
	n00=n01=n02=n10=n11=n12=n20=n21=n22=0
	n10 = 10
	n01 = 10
	n12 = 10
	n21 = 10
	counts = (n00,n01,n02,n10,n11,n12,n20,n21,n22)
	assert four_gamete_rule(counts) == 20/(2.0 * sum(counts))
	
	n10 = 10
	n01 = 10
	n12 = 10
	n21 = 20
	counts = (n00,n01,n02,n10,n11,n12,n20,n21,n22)
	assert four_gamete_rule(counts) == 20/(2.0 * sum(counts))
	
	n10 = 10
	n01 = 10
	n12 = 20
	n21 = 10
	counts = (n00,n01,n02,n10,n11,n12,n20,n21,n22)
	assert four_gamete_rule(counts) == 20/(2.0 * sum(counts))
	
	n10 = 10
	n01 = 20
	n12 = 10
	n21 = 10
	counts = (n00,n01,n02,n10,n11,n12,n20,n21,n22)
	assert four_gamete_rule(counts) == 20/(2.0 * sum(counts))
	
	n10 = 10
	n01 = 10
	n12 = 10
	n21 = 10
	counts = (n00,n01,n02,n10,n11,n12,n20,n21,n22)
	assert four_gamete_rule(counts) == 20/(2.0 * sum(counts))
	
def main():
	four_gamete_rule_test()
	args = process_input()
	vcf_file = args.vcf_file
	out_file = args.out_file
	markers = args.marker_objs
	window = args.window
	hap_thresh = args.hap_thresh
	min_maf = args.maf
	samples_file = args.samples_file
	
	samples_to_use = dict()
	samples_indicies = []
	
	#read samples from smaple id file if provided
	if samples_file is not None:
		samples_file_h = open(samples_file,"r")
		for line in samples_file_h:
			line = line.rstrip()
			samples_to_use[line] = 0
		
		samples_file_h.close()
	
	
	# 
	####
	vcf_file_handle = VariantFile(vcf_file)
	
	
	######
	# FIND OVERLAPPING SAMPLES
	######
	
	if len(samples_to_use) == 0:
		samples_indicies = range(0,len(vcf_file_handle.header.samples))
	else:
		#create list of samples_indicies
		for i in range(0,len(vcf_file_handle.header.samples)):
			sample_id = vcf_file_handle.header.samples[i]
			if sample_id in samples_to_use:
				samples_indicies.append(i)
				samples_to_use[sample_id]  = 1
		
		for sample_id, val in samples_to_use.iteritems():
			if val == 0:
				print "Sample {0} not found!".format(sample_id)
	
	print "Using {0:d} of {1:d} samples...".format(len(samples_indicies),len(vcf_file_handle.header.samples))			
				
	out_file_h = open(out_file,"w")
	out_file_h.write("MARKER,LD_BLOCK\n")
	for marker in markers:
		print "Working on marker: {0}".format(marker)
		marker_start = marker.pos - window - 1
		marker_end = marker.pos + window
		#print marker_start
		#print marker_end
		variants = [rec for rec in vcf_file_handle.fetch(marker.chrom, marker.pos - window - 1, marker.pos + window)]
		
		print "Number of variants: {0}".format(str(len(variants)))
		
		####
		# Filter by MAF for variants
		####
		print "Filtering variants..."
		variants_to_keep = []
		for variant in variants:
			samples = variant.samples
			non_missing_sample_count = 0.0
			ac = 0.0
			for i in samples_indicies:
				sample = samples[i]
				allele1 = sample["GT"][0]
				allele2 = sample["GT"][1]
				if allele1 == '.' or allele2 == '.':
					continue
				
				#COUNT ALL ALTS THE SAME
				allele1 = 1 if allele1 > 0 else 0
				allele2 = 1 if allele1 > 0 else 0 
				ac += float(allele1)
				ac += float(allele2)
				non_missing_sample_count += 1.0
			
			af = ac/(2.0 * non_missing_sample_count)
			maf = min(af, 1-af)
			
			key = Marker.format(variant.chrom,variant.pos,variant.ref,",".join(variant.alts))
			if maf >= min_maf:
				variants_to_keep.append(variant)
				#print "Keeping: {0} due to MAF of {1}".format(key,str(maf))
			#else:
			#	print "Filtering out: {0} due to MAF of {1} < {2}".format(key,str(maf), str(min_maf))
				
		print "Number of variants after filtering: {0}".format(str(len(variants_to_keep)))
		
		print "Finding marker {0}...".format(marker.key())
		###
		# Find marker in list
		###
		marker_index = 0
		for variant in variants_to_keep:
			key = Marker.format(variant.chrom,variant.pos,variant.ref,",".join(variant.alts))
			if marker.key() == key:
				print "Found {0}".format(key)
				break
			marker_index += 1
		print "Index {0}".format(marker_index)	
		
		

		
		
		####
		# Build forward haplotype block
		####
		print "Building forward haplotype block..."
		variant_of_interest = variants_to_keep[marker_index]
		forward_variants = []
		for i in range(marker_index,len(variants_to_keep) - 1):
			var_a = variants_to_keep[i]
			var_b = variants_to_keep[i + 1]
			var_gts = get_gts_from_variants(var_a,var_b,samples_indicies)
			counts = create_four_gamete_rule_counts(var_gts)
			min_freq = four_gamete_rule(counts)
			
			if min_freq > hap_thresh:
				break
			else:
				forward_variants.append(var_b)
					 
			
		backward_variants = []
		####
		# Build backward haplotype
		####
		print "Building backward haplotype block..."
		back_indices = range(1,marker_index + 1)
		back_indices.reverse()
		for i in back_indices:
			var_a = variants_to_keep[i]
			var_b = variants_to_keep[i - 1]
			var_gts = get_gts_from_variants(var_a,var_b,samples_indicies)
			counts = create_four_gamete_rule_counts(var_gts)
			min_freq = four_gamete_rule(counts)
			
			if min_freq > hap_thresh:
				break
			else:
				backward_variants.append(var_b)
	
		####
		# Merge variant keys in haplotype block
		####
		all_markers = []
		backward_variants.reverse()
		for variant in itertools.chain(backward_variants, [variant_of_interest], forward_variants):
			key = Marker.format(variant.chrom,variant.pos,variant.ref,",".join(variant.alts))
			all_markers.append(key)
		
		out_file_h.write(str(marker))
		out_file_h.write(",")
		out_file_h.write(" ".join(all_markers))
		out_file_h.write("\n")
		
	vcf_file_handle.close()
	out_file_h.close()

def create_four_gamete_rule_counts(var_gts):
	var_a_gts, var_b_gts = var_gts
	assert len(var_a_gts) == len(var_b_gts)
	n00=n01=n02=n10=n11=n12=n20=n21=n22=0
	
	for i in range(0,len(var_a_gts)):
		gt_a = var_a_gts[i]
		gt_b = var_b_gts[i]
		
		if gt_a == 0 and gt_b == 0:
			n00 += 1
		elif gt_a == 0 and gt_b == 1:
			n01 += 1
		elif gt_a == 0 and gt_b == 2:
			n02 += 1
		elif gt_a == 1 and gt_b == 0:
			n10 += 1
		elif gt_a == 1 and gt_b == 1:
			n11 += 1
		elif gt_a == 1 and gt_b == 2:
			n12 += 1
		elif gt_a == 2 and gt_b == 0:
			n20 += 1
		elif gt_a == 2 and gt_b == 1:
			n21 += 1
		elif gt_a == 2 and gt_b == 2:
			n22 += 1
		else:
			raise ValueError("Unsupported GT combination: {0} and {1}".format(gt_a,gt_b))
	
	counts = (n00,n01,n02,n10,n11,n12,n20,n21,n22)
	return counts
	
def get_gts_from_variants(var_a,var_b,samples_indicies):
	var_a_gts = []
	var_b_gts = []
	
	samples_a = var_a.samples
	for i in samples_indicies:
		gt = None
		sample = samples_a[i]
		allele1 = sample["GT"][0]
		allele2 = sample["GT"][1]
		if allele1 != '.' and allele2 != '.':
			allele1 = 1 if allele1 > 0 else 0
			allele2 = 1 if allele1 > 0 else 0
			gt = allele1 + allele2
		
		var_a_gts.append(gt)
		
	samples_b = var_b.samples
	for i in samples_indicies:
		gt = None
		sample = samples_b[i]
		allele1 = sample["GT"][0]
		allele2 = sample["GT"][1]
		if allele1 != '.' and allele2 != '.':
			allele1 = 1 if allele1 > 0 else 0
			allele2 = 1 if allele1 > 0 else 0
			gt = allele1 + allele2

		var_b_gts.append(gt)
	
	var_a_gts_no_missing = []
	var_b_gts_no_missing = []
	
	assert len(var_a_gts) == len(var_b_gts)
	
	for i in range(0,len(var_a_gts)):
		gt_a = var_a_gts[i]
		gt_b = var_b_gts[i]
		
		if gt_a is not None and gt_b is not None:
			var_a_gts_no_missing.append(gt_a)
			var_b_gts_no_missing.append(gt_b)
			
	return (var_a_gts_no_missing,var_b_gts_no_missing)
	
class Marker:

	def __init__(self,location):
		match = re.match("(chr)?([^:]+):(\d+):([ACTG]+):([ACTG]+)",location)
		if match is None:
			raise ValueError("{0} not formatted correctly. Should be chr:pos:ref:alt".format(location))
		self.chrom = match.group(2)
		self.pos = int(match.group(3))
		self.ref = match.group(4)
		self.alt = match.group(5)
	
	def key(self):
		return self.__str__()
		
	@staticmethod
	def format(chrom,pos,ref,alt):
		return "{0}:{1}:{2}:{3}".format(chrom,pos,ref,alt)	
	
	def __str__(self):
		return Marker.format(self.chrom,self.pos,self.ref,self.alt)
	
	def __repr__(self):
		return self.__str__()
		
def process_input():
	##############
	# INPUT PROCESSSING
	##############

	calling_dir = os.path.dirname(os.path.realpath(__file__))
	print "########################################################"
	print "Define variants in haplotype using Four Gamete Rule"
	print "########################################################"
	parser = argparse.ArgumentParser(description='This program defines haplotype blocks using the Four Gamete Rule')
	parser.add_argument('--vcf_file', help='the vcf file to select out of (must be bcf or tabixed)', required=True)
	parser.add_argument('--out_file', help='the location of the program output', required=True)
	parser.add_argument('--markers', help='A list of CHR:POS to define haplotype blocks for', nargs='+', type=str, required=True)
	
	parser.add_argument('--samples_file', help='Only use samples in the file provided (1 SAMPLE_ID per line)', type=str, required=False, default=None)
	parser.add_argument('--window', help='The window size in base pairs', type=int, required=False, default=50000)
	parser.add_argument('--hap_thresh', help='All haplotypes with a frequency greater than --hap_thresh will be counted', type=float, required=False, default=0.01)
	parser.add_argument('--maf', help='Only consider variants with an allele frquency greater than --maf', type=float, required=False, default=0.05)
	args = parser.parse_args()

	print

	print "########################################################"
	print "OPTIONS"
	print "########################################################"

	print

	for attr, value in args.__dict__.iteritems():
		print "{0:>25}\t\t{1!s}".format( "--{0}".format(attr), str(value))

	if len(args.markers) == 0:
		raise ValueError("Please specify at least 1 position")
	
	marker_objs = []
	for mar in args.markers:
		marker_obj = Marker(mar)
		marker_objs.append(marker_obj)
	
	args.marker_objs = marker_objs
	
	return args
	
	
	
if __name__ == '__main__':
	main()