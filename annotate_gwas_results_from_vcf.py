#!/usr/bin/python

import argparse
import re
import pysam
import os
import sys
import string
from pysam import VariantFile
import itertools
import gzip

"""
annotate_gwas_results_from_vcf.py This script will take a vcf file and gwas results and add additional columns from a 
SnpEff annotated VCF
"""

##############
# INPUT PROCESSSING
##############

calling_dir = os.path.dirname(os.path.realpath(__file__))
print "########################################################"
print "This script annotated GWAS results file from a SNPEff annotated vcf"
print "########################################################"
parser = argparse.ArgumentParser(description='This script annotated GWAS results file from a SNPEff annotated vcf')
parser.add_argument('--gwas', help='The gwas results file (tab separated)', required=True)
parser.add_argument('--vcf', help='SNPEff annotated vcf', required=True)
parser.add_argument('--out', help='The output file', required=True)
parser.add_argument('--match_by_id', help='Match the vcf entries by marker id. Make sure to specify the column number for the marker id in gwas file', type=bool, default=False, required=False)
parser.add_argument('--match_by_CHR_POS_REF_ALT', help='Match the vcf entries by CHR:POS:REF:ALT. Make sure to specify the column number for chr,pos,ref, and alt in gwas file', type=bool, default=False, required=False)
parser.add_argument('--chrom', help='Chromosome column number in gwas file', type=int, default=0, required=False)
parser.add_argument('--pos', help='Position column number in gwas file', type=int, default=1, required=False)
parser.add_argument('--marker_id', help='Marker ID column number in gwas file', type=int, default=3, required=False)
parser.add_argument('--ref', help='REF column number in gwas file', type=int, default=3, required=False)
parser.add_argument('--alt', help='ALT column number in gwas file', type=int, default=4, required=False)
parser.add_argument('--skip_lines', help='How many lines to skip in gwas results', type=int, default=0, required=False)
parser.add_argument('--info', action='append', help='Additional info fields to append to gwas results. This can be specified multiple times')

args = parser.parse_args()

print

print "########################################################"
print "OPTIONS"
print "########################################################"

print

for attr, value in args.__dict__.iteritems():
        print "{0:>15}\t\t{1!s}".format( "--{0}".format(attr), str(value))

gwas = args.gwas
vcf = args.vcf
out = args.out
match_by_id = args.match_by_id
match_by_CHR_POS_REF_ALT = args.match_by_CHR_POS_REF_ALT
chrom_ind =  args.chrom
pos_ind = args.pos
marker_id_ind = args.marker_id
ref_ind = args.ref
alt_ind = args.alt
skip_lines = args.skip_lines
infos = args.info

input_handle = None
if gwas.endswith("gz"):
	input_handle = gzip.open(gwas,"r")
else:
	input_handle = open(gwas,"r")

output_handle = None
if out.endswith("gz"):
	output_handle = gzip.open(out,"w")
else:
	output_handle = open(out,"w")

vcf_handle = None

######
# Skip lines
######
for i in range(0,skip_lines):
	input_handle.readline()

######
# Write header
######
header_cols = input_handle.readline().rstrip().split("\t")
header_cols.append("rsid")
header_cols.append("allele")
header_cols.append("gene")
header_cols.append("annotation")
header_cols.append("hgvs_c")
header_cols.append("hgvs_p")
for info in infos:
	header_cols.append(info)

output_handle.write("\t".join(header_cols))
output_handle.write("\n")


######
# Process GWAS file
######

count = 0
# create VCF handle
vcf_handle = VariantFile(vcf,"r")
	
for line in input_handle:
	count += 1
	if count % 1000 == 0:
		print "Processed {0} variants".format(count)
		
	cols = line.rstrip().split("\t")
	chrom = cols[chrom_ind]
	pos = int(cols[pos_ind])
	ref = cols[ref_ind]
	alt = cols[alt_ind]
	marker_id = cols[marker_id_ind]
	rsid = "."
	infos_out = []
	allele = "."
	gene = "."
	annotation = "."
	hgvs_c = "."
	hgvs_p = "."
	
	#####
	# Get annotation
	#####
	for rec in vcf_handle.fetch(chrom,pos - 1, pos + 1):
		####
		# Match by position
		####
		if rec.chrom == chrom and rec.pos == pos:
			
			if "ANN" in rec.info:
				ann_field = rec.info["ANN"]
					#anns = ann_field.split(",")
				ann = ann_field[0]
				ann_cols = ann.split("|")
				allele = ann_cols[0]
				gene = ann_cols[3]
				annotation = ann_cols[1]
				hgvs_c = ann_cols[9]
				hgvs_p = ann_cols[10]
				
			rsid = rec.id

			for info in infos:
				if info in rec.info:
					infos_out.append(rec.info[info])
				else:
					infos_out.append(".")
			
			
			# leave the loop
			break

	#####
	# Write out annotation
	#####
	
	# make surre info is full
	if len(infos_out) == 0:
		for info in infos:
			infos_out.append(".")
	
	###
	# Add annotation to cols and write to file
	###
	cols.append(str(rsid))
	cols.append(str(allele))
	cols.append(str(gene))
	cols.append(str(annotation))
	cols.append(str(hgvs_c))
	cols.append(str(hgvs_p))
	for info in infos_out:
		if type(info) is tuple or type(info) is list:
			info_str = [str(s) for s in info]
			cols.append(",".join(info_str))
		else:
			cols.append(str(info))
	
	output_handle.write("\t".join(cols))
	output_handle.write("\n")
	
######
# Clean up
######
if input_handle is not None:
	input_handle.close()

if output_handle is not None:
	output_handle.close()

if vcf_handle is not None:
	vcf_handle.close()

print "Complete!"