#!/usr/bin/python

import argparse
import re
import pysam
import pickle
import os
import string
from pysam import VariantFile
import itertools
#import gzip
	
def main():
	args = process_input()
	index_dbsnp = args.index_dbsnp
	vcf_file = args.vcf_file
	out_file = args.out_file
	marker_ids = args.marker_ids
	dbsnp_index = args.dbsnp_index
	
	if index_dbsnp:
		index_dbsnp_func(vcf_file,out_file)
	else:
		create_haplotypes(marker_ids,vcf_file,out_file,dbsnp_index)


def process_input():
	##############
	# INPUT PROCESSSING
	##############

	calling_dir = os.path.dirname(os.path.realpath(__file__))
	print "########################################################"
	print "Extract haplotypes from vcf"
	print "########################################################"
	parser = argparse.ArgumentParser(description='This program creates a haplotype for each subject from minimac results')
	parser.add_argument('--out_file', help='the output file', required=True)
	parser.add_argument('--vcf_file', help='the vcf file to select out of (must be bcf or tabixed)', required=True)

	parser.add_argument('--dbsnp_index', help='dbsnp index database', required=False, default=None)
	parser.add_argument('--marker_ids', help='the ids of create haploptyes for', nargs='+', type=str, required=False)
	parser.add_argument('--index_dbsnp', help='flag to index the dbsnp database', required=False, action='store_true')
	args = parser.parse_args()

	print

	print "########################################################"
	print "OPTIONS"
	print "########################################################"

	print

	for attr, value in args.__dict__.iteritems():
		print "{0:>25}\t\t{1!s}".format( "--{0}".format(attr), str(value))
		
	return args

class DbSnpEntry:
	
	def __init__(self,rs,chrom,pos):
		self.rs_int = long(rs.lstrip("rs")) 
		self.chrom = chrom
		self.pos = pos
	
	def rs(self):
		return "rs{0}".format(self.rs_int)

	def __str__(self):
		return "{0}, {1}:{2}".format(self.rs(), self.chrom,self.pos)
		
def index_dbsnp_func(dbsnp_vcf,out_file):
	print "########################################################"
	print "INDEXING DBSNP"
	print "########################################################"
	vcf_in = VariantFile(dbsnp_vcf)
	count = 0
	limit = 1000000
	rs_dict = dict()
	for rec in vcf_in.fetch():
		count += 1
		
		if count > limit:
			break;
			
		rs = rec.id
		chrom = rec.chrom
		pos = rec.pos
		
		entry = DbSnpEntry(rs,chrom,pos)
		
		if count % 100000 == 0:
			print "LINE#{0} {1}".format(count,entry)
					
		rs_dict[entry.rs_int] = entry
	
	print "WRITING OUT dbsnp index... to {0}".format(out_file)
	pickle.dump(rs_dict, open(out_file,'wb') )

class ChrPos:
	
	def __init__(self,chrom,pos):
		self.chrom = chrom
		self.pos = pos
	
	def __str__(self):
		return "{0}:{1}".format(self.chrom,self.pos)

class SampleHaplotype:
	
	def __init__(self,ident):
		self.id = ident
		self.haplotype_1 = []
		self.haplotype_2 = []
	
	def add_allele_to_haplotype_1(self,allele):
		self.haplotype_1.append(allele)
	
	def add_allele_to_haplotype_2(self,allele):
		self.haplotype_2.append(allele)
				
	def get_haplotype_1(self):
		return "".join(self.haplotype_1)
		
	def get_haplotype_2(self):
		return "".join(self.haplotype_2)
			
	def get_haplotypes(self):
		return [self.get_haplotype_1(), self.get_haplotype_2()]
	
	def dosage(self,haplotype):
		count = 0
		if haplotype == self.get_haplotype_1():
			count += 1
			
		if haplotype == self.get_haplotype_2():
			count += 1
		
		return count
		
		
	def __str__(self):
		return "{id}: {hap1}/{hap2}".format(id=self.id, hap1=self.get_haplotype_1(), hap2=self.get_haplotype_2())

def create_haplotypes(marker_ids,vcf_file,out_file,dbsnp_index=None):
	print "########################################################"
	print "CREATING HAPLOTYPES"
	print "########################################################"
	
	vcf_file_handle = VariantFile(vcf_file)
	
	dbsnp_dict = None
	if dbsnp_index is not None:
		print "Reading dbsnp..."
		dbsnp_dict = pickle.load(open(dbsnp_index,'rb'))
		print "Read dbsnp"
	
	chr_pos_vals = []
	
	for marker in marker_ids:
		if marker.startswith("rs"):
			raise NotImplementedError("Only chr:pos supported currently")
		else:
			chrom,pos = marker.split(":")
			pos = int(pos)
			chr_pos = ChrPos(chrom,pos)
			chr_pos_vals.append(chr_pos)
	
	###
	# Make sure all markers are on the same chr
	###
	chroms = dict()
	for chr_pos in chr_pos_vals:
		chroms[chr_pos.chrom] = 1
	
	if len(chroms.keys()) > 1:
		raise ValueError("Markers are not all on same chromosome")
	
	print "Input appears valid..."
	print "Reading genotype data at each site.."
	gts_for_markers = dict()
	sample_ids = [x for x in vcf_file_handle.header.samples ]
	for chr_pos in chr_pos_vals:
		for rec in vcf_file_handle.fetch(chr_pos.chrom, chr_pos.pos - 1, chr_pos.pos):
			marker_key = "{0}:{1}:{2}:{3}".format(rec.chrom,rec.pos,rec.ref,"_".join(rec.alts))
			print "Reading Marker: {0}".format(marker_key)
			if rec.pos == chr_pos.pos:
				marker_dict = dict()
				samples = rec.samples
				for i in range(0,len(samples)):
					sample = samples[i]
					if sample.phased:
						marker_dict[sample.name] = sample["GT"]
					else:
						marker_dict[sample.name] = [".", "."]
			gts_for_markers[marker_key] = marker_dict
	vcf_file_handle.close()
	marker_keys = gts_for_markers.keys()
	
	###
	# Discover haplotypes
	###
	sample_haplotypes = dict()
	haplotypes = dict()
	print "Markers:"
	print marker_keys
	for sample_id in sample_ids:
		sample_haplotype = SampleHaplotype(sample_id)
		for maker_id in marker_keys:
			gt_dict = gts_for_markers[maker_id]
			if sample_id in gt_dict.keys():
				allele_1 = str(gt_dict[sample_id][0])
				allele_2 = str(gt_dict[sample_id][1])
				sample_haplotype.add_allele_to_haplotype_1(allele_1)
				sample_haplotype.add_allele_to_haplotype_2(allele_2)
			
		sample_haplotypes[sample_haplotype.id] = sample_haplotype
	
	for sample_id in sample_ids:
		sample_haplotype = sample_haplotypes[sample_id]
		if sample_haplotype.get_haplotype_1() in haplotypes.keys():
			haplotypes[sample_haplotype.get_haplotype_1()] += 1
		else:
			haplotypes[sample_haplotype.get_haplotype_1()] = 1
			
		if sample_haplotype.get_haplotype_2() in haplotypes:
			haplotypes[sample_haplotype.get_haplotype_2()] += 1
		else:
			haplotypes[sample_haplotype.get_haplotype_2()] = 1
	
	print "Haplotypes discovered:"
	haplotype_tot = float(2 * len(sample_ids))
	for haplotype,count in haplotypes.iteritems():
		print "{0}: {1}".format(haplotype, count / haplotype_tot)
	
	print "Writing out haplotype table..."
	out_file_handle = open(out_file,"w")
	hap_pos = 0
	for maker_id in marker_keys:
		hap_pos += 1
		chrom,pos,ref,alts = maker_id.split(":")
		out_file_handle.write("#{0}:{1}\tHAP_POS:{2}\t".format(chrom,pos,hap_pos))
		out_file_handle.write("0 = {0}, ".format(ref))
		out_file_handle.write(" 1 = {0}".format(alts))
		out_file_handle.write("\n")
	
	hap_keys = haplotypes.keys()
	header = itertools.chain(["SAMPLE_ID"],hap_keys)
	
	
	out_file_handle.write("\t".join(header))
	out_file_handle.write("\n")
	
	for sample_id in sample_ids:
		cols = [sample_id]
		sample_haplotype = sample_haplotypes[sample_id]
		for haplotype in hap_keys:
			cols.append( str( sample_haplotype.dosage(haplotype)))
		
		out_file_handle.write("\t".join(cols))
		out_file_handle.write("\n")
		
	out_file_handle.close()
	print "Complete!"
	
if __name__ == '__main__':
	main()
		