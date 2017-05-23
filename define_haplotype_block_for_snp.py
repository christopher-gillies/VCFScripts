#!/usr/bin/python

import argparse
import re
import os
import string
from pysam import VariantFile
import itertools
	
def main():
	args = process_input()
	vcf_file = args.vcf_file
	out_file = args.out_file
	markers = args.marker_objs
	window = args.window
	hap_freq = args.hap_freq
	maf = args.maf
	
	####
	#
	####
	vcf_file_handle = VariantFile(vcf_file)
	for marker in markers:
		variants = [rec for rec in vcf_file_handle.fetch(marker.chrom, marker.pos - 1 - window, marker.pos + window)]
		
		#Find index of marker of interest in list
		variant_markers = []
		for variant in variants:
			raise NotImplementedError()
		###
		# Find marker in list
		###
	
class Marker:

	def __init__(self,location):
		match = re.match("(chr)([^:]+):(\d+):([ACTG]+):([ACTG]+)")
		if match is None:
			raise ValueError("{0} not formatted correctly. Should be chr:pos:ref:alt".format(location))
		self.chrom = match.group(1)
		self.pos = int(match.group(2))
		self.ref = match.group(3)
		self.ref = match.group(4)
		
	def __str__(self):
		return "{0}:{1}:{2}:{3}".format(self.chrom,self.pos,self.ref,self.alt)	
		
	
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
	
	parser.add_argument('--window', help='The window size in base pairs', nargs='+', type=int, required=False, default=50000)
	parser.add_argument('--hap_freq', help='All haplotypes with a frequency greater than --hap_freq will be counted', nargs='+', type=float, required=False, default=0.01)
	parser.add_argument('--maf', help='Only consider variants with an allele frquency greater than --maf', nargs='+', type=float, required=False, default=0.05)
	args = parser.parse_args()

	print

	print "########################################################"
	print "OPTIONS"
	print "########################################################"

	print

	for attr, value in args.__dict__.iteritems():
		print "{0:>25}\t\t{1!s}".format( "--{0}".format(attr), str(value))

	if len(positions) == 0:
		raise ValueError("Please specify at least 1 position")
	
	marker_objs = []
	for mar in args.markers:
		marker_obj = Marker(pos)
		marker_objs.append(marker_obj)
	
	args.marker_objs = marker_objs
	
	return args
	
	
	
if __name__ == '__main__':
	main()