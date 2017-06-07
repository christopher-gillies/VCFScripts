#!/usr/bin/python

import argparse
import re
import os
import pysam
#from subprocess import Popen, PIPE
from pysam import VariantFile


def process_input():
	##############
	# INPUT PROCESSSING
	##############
	
	calling_dir = os.path.dirname(os.path.realpath(__file__))
	print "########################################################"
	print "Filters imputed vcf by r2"
	print "########################################################"
	parser = argparse.ArgumentParser(description='Convert Illumina standard reports to lgen format')
	parser.add_argument('--vcfs', help='a file containing vcfs. each line should contain chr [tab] vcf', required=True)
	parser.add_argument('--out_prefix', help='the out prefix for the vcf files', required=True)

	parser.add_argument('--min_r2', help='the minimum r2 to keep', type=float, required=False, default = 0.7)
	parser.add_argument('--min_maf', help='the minimum maf to keep', type=float, required=False, default = 0.01)
	parser.add_argument('--r2_field_name', help='default field name is R2', type=str, required=False, default = "R2")
	parser.add_argument('--maf_field_name', help='default field name is MAF', type=str, required=False, default = "MAF")
	parser.add_argument('--new_ids', help='A two column file containing OLD_ID and NEW_ID and a new line. If an id is not present, it will not be relabelled', type=str, required=False, default = None)
	args = parser.parse_args()


	print

	print "########################################################"
	print "OPTIONS"
	print "########################################################"

	print

	for attr, value in args.__dict__.iteritems():
		print "{0:>25}\t\t{1!s}".format( "--{0}".format(attr), str(value))
	
	print "Reading vcfs file list: {0}".format(args.vcfs)
	
	chrom_vcf = dict()
	with open(args.vcfs,"r") as vcfs:
		for line in vcfs:
			chrom,vcf = line.rstrip().split("\t")
			chrom_vcf[chrom] = vcf
	
	args.chrom_vcf = chrom_vcf
	
	return args
	

def main():
	args = process_input()
	
	chrom_vcf = args.chrom_vcf
	min_r2 = args.min_r2
	min_maf = args.min_maf
	out_prefix = args.out_prefix
	r2_field_name = args.r2_field_name
	maf_field_name = args.maf_field_name
	new_ids = args.new_ids
	
	####
	# Read new ids in dictionary
	####
	
	new_ids_dict = dict()
	if new_ids is not None:
		with open(new_ids,"r") as f:
			for line in f:
				old_id,new_id = line.rstrip().split("\t")
				new_ids_dict[old_id] = new_id
		print "Ids {0} ids to remap".format(len(new_ids_dict))		
				
	out_vcf_list = "{0}.vcf_list.tsv".format(out_prefix)
	out_vcf_list_handle = open(out_vcf_list,"w")
	
	
	
	for chrom,vcf in chrom_vcf.iteritems():
		chrom_match = re.match("(chr)?(.+)", chrom)
		if chrom_match is not None:
			chrom = chrom_match.group(2)
		else:
			raise ValueError("Chomosome name {0} not formatted correctly!".format(chrom))
		
		
		out_vcf_name = "{0}.chr{1}.vcf".format(out_prefix,chrom)
		out_vcf_name_gz = "{0}.chr{1}.vcf.gz".format(out_prefix,chrom)
		out_vcf_name_gz_tbi = "{0}.chr{1}.vcf.gz.tbi".format(out_prefix,chrom)
		
		print "Processing chr{0} {1}...".format(chrom,vcf)
		in_vcf_handle = VariantFile(vcf)
		pass_filter = in_vcf_handle.header.filters["PASS"]
		
		out_vcf_list_handle.write("{0}\t{1}".format(chrom,out_vcf_name_gz))
		out_vcf_list_handle.write("\n")
		
		####
		# It appears that writing to a BCF is the only method that works in this version of pysam
		####
		
		#'wb' for BCF
		#
		#out_vcf_handle = VariantFile(out_vcf_name,'wb',header=in_vcf_handle.header)
		#out_vcf_handle = pysam.libcbgzf.BGZFile(out_vcf_name,"wb")
		#out_vcf_handle.write(str(in_vcf_handle.header))
		
		#cmd = "bgzip -c > {0}".format(out_vcf_name)
		#print cmd
		
		
		out_vcf_handle = open(out_vcf_name,"w")
		
		print "Relabeling and writing header..."
		relabeled_ids = 0
		old_header_lines = str(in_vcf_handle.header).split("\n")
		for line in old_header_lines:
			
			if line == "":
				continue
			
			if re.match("^#CHROM.+",line):
				cols = line.split("\t")
				for i in range(9,len(cols)):
					if cols[i] in new_ids_dict:
						relabeled_ids += 1
						cols[i] = new_ids_dict[cols[i]]
				#merge new columns
				new_line = "\t".join(cols)
				out_vcf_handle.write(new_line)
			else:
				out_vcf_handle.write(line)
			
			#write new line
			out_vcf_handle.write("\n")
		
		print "Relabeled {0} ids".format(relabeled_ids)
		
		rec_count = 0

		for rec in in_vcf_handle:
			rec_count += 1
			if rec_count % 50000 == 0:
				print "Line: {0:d} {1}:{2:d}".format(rec_count,rec.chrom,rec.pos)
			r2 = rec.info[r2_field_name]
			maf = rec.info[maf_field_name]
			if r2 > min_r2 and maf > min_maf:
				#clear filters
				rec.filter.clear()
				#set filter to be pass
				rec.filter.add("PASS")
				#new lines are already there
				out_vcf_handle.write(str(rec))
				
		#print "Running bgzip on "
		##execute bgzip
		#bgz_handle = Popen(["bgzip", out_vcf_name])
		#bgz_handle.wait()
		
		in_vcf_handle.close()
		out_vcf_handle.close()
		
		print "Writing tabix index for {0}...".format(out_vcf_name, preset="vcf")
		#seems to only compress files
		pysam.tabix_index(out_vcf_name,preset="vcf")
		
		if not os.path.isfile(out_vcf_name_gz_tbi):
			pysam.tabix_index(out_vcf_name_gz,preset="vcf")
		
		if os.path.isfile(out_vcf_name):
			os.remove(out_vcf_name)
			
	out_vcf_list_handle.close()
	print "Finished writing {0}".format(out_vcf_list)
	print "Complete!"

if __name__ == '__main__':
	main()
