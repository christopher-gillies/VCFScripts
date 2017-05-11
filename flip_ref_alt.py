#!/usr/bin/python

import argparse
import re
import gzip
import sys
import string
parser = argparse.ArgumentParser(description='Take results from checkVCF.py reference flips and flip ref and alt')
parser.add_argument('--vcf', help='vcf file (bgzip or vcf)', required=True)
parser.add_argument('--ref_flips', help='ref flip file from checkVCF.py', required=True)

args = parser.parse_args()

vcf = args.vcf
ref_flips =  args.ref_flips


sys.stderr.write("Input VCF Format: ")
gzipped = False
if re.search("vcf.gz$", vcf) != None:
	sys.stderr.write("GZIP\n")
	gzipped = True
elif re.search("vcf$", vcf) != None:
	sys.stderr.write("STD_VCF\n")
else:
	sys.stderr.write("Not a VCF file\n")
	quit()
	

####
# VariantKey class takes input from checkVCF.py
####

class VariantKey:
	
	def __init__(self,variant):
		match = re.match("([^:]+)[:](\d+)[:]([ACTGactg])[-]([ACTGactg])[\/]([ACTGactg])", variant)
		self.valid = False
		if match != None:
			self.chr = match.group(1)
			self.pos = match.group(2)
			self.new_ref = match.group(3)
			self.old_ref = match.group(4)
			self.old_alt = match.group(5)
			self.new_alt = self.old_ref
			if self.new_ref == self.old_alt:
				self.valid = True
			self.old_key = '{0}:{1}:{2}:{3}'.format(self.chr,self.pos,self.old_ref,self.old_alt)
			self.new_key = '{0}:{1}:{2}:{3}'.format(self.chr,self.pos,self.new_ref,self.new_alt)
	
	def __str__(self):
		valid_text = "NOT VALID"
		if self.valid:
			valid_text = "VALID"
			
		return '{0} --> {1} [{2}]'.format(self.old_key, self.new_key, valid_text )


################
# Read flips file and store in dictionary
################
var_keys = dict()				
with open(ref_flips,'r') as ref_flips_handle:
	for line in ref_flips_handle:
		line = line.rstrip()
		(info,variant) = re.split("\s",line)
		variant_key = VariantKey(variant)
		if variant_key.valid:
			var_keys[variant_key.old_key] = variant_key

sys.stderr.write("Variants to swap: ")
sys.stderr.write(str(len(var_keys)))
sys.stderr.write("\n")

vcf_file_handle = None

if gzipped:
	vcf_file_handle = gzip.open(vcf, 'rb')
else:
	vcf_file_handle = open(vcf, 'r')

count = 0


###########
# Read VCF
##########

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ID1
swap_count = 0
for line in vcf_file_handle:
	count += 1

	line = line.rstrip()
	#Just print out header lines
	if re.search("^#", line) != None:
		print line
		continue
	
	cols = re.split("\t",line)
	
	#Just print multialleleic
	if re.search(",", cols[4]) != None:
		print line
		continue
		
	var_key = '{0}:{1}:{2}:{3}'.format(cols[0],cols[1],cols[3],cols[4])
	if count % 10000 == 0:
		sys.stderr.write("Line: {0} Variant: {1}\n".format(count,var_key))
		sys.stderr.write("REF/ALTs swapped: {0}\n".format(swap_count))
	if var_key in var_keys:
		swap_count += 1
		# FLIP REF AND ALT
		variant_key_obj = var_keys[var_key]
		cols[3] = variant_key_obj.new_ref
		cols[4] = variant_key_obj.new_alt
		for i in range(9,len(cols)):
			col = cols[i]
			fields = re.split("[:]", col)
			gts = fields[0]
			#convert 0 ==> 1 and 1 ==> 0
			gts_reformat = gts.translate(string.maketrans("01", "10"))
			fields[0] = gts_reformat
			#reassign gt column with new fields
			
			cols[i] = ":".join(fields)
			#reassign line with modified columns
			old_line = line
			line = "\t".join(cols)
			
			if len(old_line) != len(line):
				sys.stderr.write("Error in processing line:\n")
				sys.stderr.write(old_line)
				sys.stderr.write("\n")
				sys.exit()
			
	print line
	

vcf_file_handle.close()
sys.stderr.write("Complete!\n")
sys.stderr.write("REF/ALTs swapped: {0}\n".format(swap_count))