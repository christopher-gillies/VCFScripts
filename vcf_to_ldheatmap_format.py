#!/usr/bin/python

import argparse
import re
import pysam
import pickle
import os
import string
from pysam import VariantFile
import itertools
from marker_util import Marker

"""
vcf_to_ldheatmap_format.py outputs vcf snps in a format acceptable for LDheatmap R package
https://cran.r-project.org/web/packages/LDheatmap/vignettes/LDheatmap.pdf
"""


def process_input():
	##############
	# INPUT PROCESSSING
	##############
	
	calling_dir = os.path.dirname(os.path.realpath(__file__))
	print "########################################################"
	print "Extract haplotypes from vcf"
	print "########################################################"
	parser = argparse.ArgumentParser(description='This program outputs two files for LDheatmap a SNP genotype file and their physical distances')
	parser.add_argument('--out', help='the output prefix', required=True)
	parser.add_argument('--vcf_file', help='the vcf file to select out of (must be bcf or tabixed)', required=True)
	parser.add_argument('--marker_ids', help='chr:pos:ref:alt separated by a space', nargs='+', type=str, required=False)
	args = parser.parse_args()
	
	
	print "########################################################"
	print "OPTIONS"
	print "########################################################"

	print

	for attr, value in args.__dict__.iteritems():
		print "{0:>25}\t\t{1!s}".format( "--{0}".format(attr), str(value))

	if len(args.marker_ids) == 0:
		raise ValueError("Please specify at least 1 marker")
		
	marker_objs = []
	for mar in args.marker_ids:
		marker_obj = Marker(mar)
		marker_objs.append(marker_obj)
	
	args.marker_objs = marker_objs
	
	return args

def main():
	
	args = process_input()
	vcf_file = args.vcf_file
	out_prefix = args.out
	
	out_snps = "{0}.ldheatmap.snps".format(out_prefix)
	# SNP in col
	# IND in row
	
	out_pos = "{0}.ldheatmap.dists".format(out_prefix)
	# a position for each SNP
	markers = args.marker_objs
	
	vcf_file_handle = VariantFile(vcf_file)
	sample_ids = vcf_file_handle.header.samples
	
	### SMAPLE_ID --> MARKER_ID --> VAL
	sample_vals_by_marker = dict()
	
	# create samples in dictionary
	
	for sample_id in sample_ids:
		sample_vals_by_marker[sample_id] = dict()
		
	markers_keys = []
	positions = []
	
	for marker in markers:
		print "Working on marker: {0}".format(marker)
		marker_start = marker.pos - 1
		marker_end = marker.pos + 1
		variants = [rec for rec in vcf_file_handle.fetch(marker.chrom, marker_start, marker_end)]
		variant_match = None
		variant_match_key = None
		# get matching variant
		for variant in variants:
			variant_match = None
			key = Marker.format(variant.chrom,variant.pos,variant.ref,",".join(variant.alts))
			if key == marker.key():
				variant_match = variant
				variant_match_key = key
		
		if variant_match is not None:
			markers_keys.append(marker.key())
			positions.append(str(marker.pos))
			
			
			print "Variant match found: {0}".format(variant_match_key)
			for i in range(0,len(sample_ids)):
				sample_id = sample_ids[i]
				sample = variant_match.samples[i]
				value = "NA"
				if sample.phased and sample["GT"][0] != "." and sample["GT"][1] != ".":
					allele1 = variant_match.alleles[ sample["GT"][0] ]
					allele2 = variant_match.alleles[ sample["GT"][1] ]
					value = "/".join((allele1,allele2))
				
				sample_vals_by_marker[sample_id][key] = value
	
	
	print "Writing out output files..."
	
	###
	# Write out position file
	###
	out_pos_h = open(out_pos,"w")
	out_pos_h.write("\t".join(positions))
	out_pos_h.write("\n")
	out_pos_h.close()
	
	###
	# Write out snp file
	###
	out_snps_h = open(out_snps,"w")
	
	### header
	out_snps_h.write("SAMPLE_ID\t")
	out_snps_h.write("\t".join(markers_keys))
	out_snps_h.write("\n")
	
	#####
	# Write out data
	#####
	for sample_id in sample_ids:
		out_snps_h.write(sample_id)
		out_snps_h.write("\t")
		vals = []
		sample_markers = sample_vals_by_marker[sample_id]
		for marker_key in markers_keys:
			if marker_key in sample_markers:
				val = sample_markers[marker_key]
				vals.append(val)
			else:
				print "KEY {0} not found for {1}".format(marker_key,sample_id)
				vals.append("NA")
		out_snps_h.write("\t".join(vals))
		out_snps_h.write("\n")
		
	out_snps_h.close()
	
	
	
if __name__ == '__main__':
	main()




