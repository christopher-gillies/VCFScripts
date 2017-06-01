#!/usr/bin/python

import argparse
import re
import os

class Sample:
	
	def __init__(self,id,file):
		self.id = id
		self.file = file
	
	

def process_input():
	##############
	# INPUT PROCESSSING
	##############
	
	calling_dir = os.path.dirname(os.path.realpath(__file__))
	print "########################################################"
	print "Take a list of samples and Illumina standard reports and convert"
	print "to lgen format which can then be converted to a bed file"
	print "########################################################"
	parser = argparse.ArgumentParser(description='Convert Illumina standard reports to lgen format')
	parser.add_argument('--out_prefix', help='the output prefix', required=True)
	parser.add_argument('--sample_file', help='A tab separated file. Column 1 is sample ID and column 2 is the file path', type=str, required=True)
	parser.add_argument('--manifest_csv', help='Illumina manifest file in column separated format', type=str, required=True)
	
	parser.add_argument('--skip_indels', help='Skip indel variants in manifest', action='store_true', required=False)
	parser.add_argument('--indel_regex', help='Indel regex definition', type=str, default="[DI]", required=False)
	parser.add_argument('--snp_col_name_manifest', help='SNP Name in manifest', type=str, default="Name", required=False)
	parser.add_argument('--snp_change_col_name_manifest', help='Typically formated [A/G]', type=str, default="SNP", required=False)
	parser.add_argument('--chrom_col_name_manifest', help='Chrom column name in manifest', type=str, default="Chr", required=False)
	parser.add_argument('--pos_col_name_manifest', help='Position column name in manifest', type=str, default="MapInfo", required=False)
	parser.add_argument('--snp_col_name_report', help='SNP Name in report', type=str, default="SNP Name", required=False)
	parser.add_argument('--allele_1_col_name_report', help='allele 1 column name in report', type=str, default="Allele1 - Forward", required=False)
	parser.add_argument('--allele_2_col_name_report', help='allele 2 column name in report', type=str, default="Allele2 - Forward", required=False)

	
	parser.add_argument('--sep_report', help='delimiter default is a tab', type=str, default="\t", required=False)
	args = parser.parse_args()


	print

	print "########################################################"
	print "OPTIONS"
	print "########################################################"

	print

	for attr, value in args.__dict__.iteritems():
		print "{0:>25}\t\t{1!s}".format( "--{0}".format(attr), str(value))
	
	sample_file_h = open(args.sample_file,"r")
	samples = []
	for line in sample_file_h:
		cols = line.rstrip().split("\t")
		sample = Sample(cols[0],cols[1])
		samples.append(sample)
		
	sample_file_h.close()
	
	args.samples = samples
	
	return args
	
def main():
	args = process_input()
	
	skip_indels = args.skip_indels
	indel_regex = args.indel_regex
	
	samples = args.samples
	manifest_csv = args.manifest_csv
	snp_col_name_manifest = args.snp_col_name_manifest
	snp_change_col_name_manifest = args.snp_change_col_name_manifest
	chrom_col_name_manifest = args.chrom_col_name_manifest
	pos_col_name_manifest = args.pos_col_name_manifest
	
	snp_col_name_report = args.snp_col_name_report
	allele_1_col_name_report = args.allele_1_col_name_report
	allele_2_col_name_report = args.allele_2_col_name_report
	
	out_prefix = args.out_prefix
	
	map_out = "{0}.map".format(out_prefix)

	sep_report = args.sep_report

	
	"""
	MAP COLS:
		 chromosome (1-22, X, Y or 0 if unplaced)
	     rs# or snp identifier
	     Genetic distance (morgans)
	     Base-pair position (bp units)
	"""
	lgen_out = "{0}.lgen".format(out_prefix)
	"""
	LGEN COLS:
	
		 family ID
	     individual ID
	     snp ID
	     allele 1 of this genotype
	     allele 2 of this genotype
	"""
	fam_out = "{0}.fam".format(out_prefix)
	
	"""
	FAM COLS:
	first six columns of PED
	
		 Family ID
	     Individual ID
	     Paternal ID
	     Maternal ID
	     Sex (1=male; 2=female; other=unknown)
	     Phenotype
	"""
	
	print "Checking that all files exists"
	
	if not os.path.isfile(manifest_csv):
		print "Manifest {0} does not exists!".format(manifest_csv)
	
	for sample in samples:
		if not os.path.isfile(sample.file):
			print "Sample {0} file {1} does not exists!".format(sample.id,sample.file)
	
	print "All files appear to exists."
	
	print "Writing out FAM file: {0}".format(fam_out)
	
	fam_out_h = open(fam_out,"w")
	for sample in samples:
		line_out = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(sample.id,sample.id,0,0,0,0)
		fam_out_h.write(line_out)
	fam_out_h.close()
	
	print "Finished writing FAM file: {0}".format(fam_out)
	
	print 
	print "Processing manifest file {0}...".format(manifest_csv)
	print 
	
	
	manifest_h = open(manifest_csv,"r")
	map_out_h = open(map_out,"w")
	# Disover header
	
	
	snp_col_number_manifest = None
	snp_change_col_number_manifest = None
	chrom_col_number_manifest = None
	pos_col_number_manifest = None
	manifest_header_col_number = None
	
	###
	# SNPs to keep from manifest
	###
	snps_to_keep = dict()
	
	manifest_line_num = 0
	
	for line in manifest_h:
		manifest_line_num += 1
		cols = line.rstrip().split(",")
		if (snp_col_name_manifest in cols and 
			snp_change_col_name_manifest in cols and chrom_col_name_manifest in cols and
			pos_col_name_manifest in cols):
			print "Header for manifest discovered..."
			snp_col_number_manifest = cols.index(snp_col_name_manifest)
			snp_change_col_number_manifest = cols.index(snp_change_col_name_manifest)
			chrom_col_number_manifest = cols.index(chrom_col_name_manifest)
			pos_col_number_manifest = cols.index(pos_col_name_manifest)
			manifest_header_col_number = len(cols)
			break
	
	if snp_col_number_manifest is None:
		raise ValueError("Header not discovered in manifest file please check manifest column names...")
	
	for line in manifest_h:
		manifest_line_num += 1
		cols = line.rstrip().split(",")
		
		if len(cols) != manifest_header_col_number:
			print "Number of columns no longer match header at line {0:d}".format(manifest_line_num)
			break
		
		snp_id = cols[snp_col_number_manifest]
		chrom = cols[chrom_col_number_manifest]
		pos = cols[pos_col_number_manifest]
		snp_change = cols[snp_change_col_number_manifest]
		
		if manifest_line_num % 50000 == 0:
			print "At SNP: {0} and line: {1}".format(snp_id,manifest_line_num)
			
		match = re.search(indel_regex,snp_change)
		
		if skip_indels and match is not None:
			continue
		else:
			# Write to marker file
			map_out_line = "{chrom}\t{snp_id}\t{morgan}\t{pos}\n".format(chrom=chrom,snp_id=snp_id,morgan=0,pos=pos)
			map_out_h.write(map_out_line)
			snps_to_keep[snp_id] = None
	
	manifest_h.close()
	map_out_h.close()
	print "Finished processing manifest file"
	print "Finished writing: {0}".format(map_out)
	print "Found {0} variants in manifest".format(len(snps_to_keep))
	
	
	print
	print "Processing {0:d} samples...".format(len(samples))
	print
	
	marker_lines = []
	header = None
	
	lgen_out_h = open(lgen_out,"w")
	
	
	sample_count = 0
	snp_col_number_report = None
	allele_1_col_number_report = None
	allele_2_col_number_report = None
	
	for sample in samples:
		sample_count += 1
		
		
		print "Processing {0}".format(sample.id)
		print "Processing file: {0}".format(sample.file)
		file_handle = open(sample.file,"r")
		
		#Get to first data line
		lines_skipped = 0
		for line in file_handle:
			line = line.rstrip()
			lines_skipped += 1
			#Match header line
			cols = line.split(sep_report)
			if (snp_col_name_report in cols and
				allele_1_col_name_report in cols and
				allele_2_col_name_report in cols):
				
				snp_col_number_report = cols.index(snp_col_name_report)
				allele_1_col_number_report = cols.index(allele_1_col_name_report)
				allele_2_col_number_report = cols.index(allele_2_col_name_report)
				print "Header:"
				print line
				break
		
		if snp_col_number_report is None:
			raise ValueError("Header not discovered in report file please check report column names...")
		
		print "Lines skipped for header: {0:d}".format(lines_skipped)
		
		# read data
		line_count = 0
		for line in file_handle:
			#print line
			cols = line.rstrip().split(sep_report)
			snp_id = cols[snp_col_number_report]
			allele_1 = cols[allele_1_col_number_report]
			allele_2 = cols[allele_2_col_number_report]
			
			if line_count % 50000 == 0:
				print "At SNP: {0} and line: {1} for sample: {2}".format(snp_id,line_count,sample.id)
			
			#skip
			if snp_id not in snps_to_keep:
				continue
			
			if allele_1 == "-":
				allele_1 = "0"
			
			if allele_2 == "-":
				allele_2 = "0"
			
			lgen_line = "{fam_id}\t{ind_id}\t{snp_id}\t{allele_1}\t{allele_2}\n".format(fam_id=sample.id, ind_id=sample.id, snp_id=snp_id,allele_1=allele_1,allele_2=allele_2)
			lgen_out_h.write(lgen_line)
			line_count += 1
		
		print "Lines {0} processed for sample {1}".format(line_count,sample.id)
			
		file_handle.close()
	
	lgen_out_h.close()
	print "Finished writing: {0}".format(lgen_out)
	
	print 
	print "OUTPUT FILES"
	print "FAM: {0}".format(fam_out)
	print "LGEN: {0}".format(lgen_out)
	print "MAP: {0}".format(map_out)
	print "Run: plink --lfile {0} --make-bed --out {0}".format(out_prefix)
	
if __name__ == '__main__':
	main()
