#!/usr/bin/python

import argparse
import re
import gzip
import sys
import string
import os
from make import MakeEntry
from make import MakeFile

##############
# INPUT PROCESSSING
##############

calling_dir = os.path.dirname(os.path.realpath(__file__))
print "########################################################"
print "Create imputation pipeline"
print "########################################################"
parser = argparse.ArgumentParser(description='This program creates and imputation pipeline for input vcfs. Imputation is performed using minimac3 and phasing is done using Eagle')
parser.add_argument('--outdir', help='the output directory', required=True)
parser.add_argument('--vcf_files', help='the list of vcf files [chr\tvcf]', required=True)
parser.add_argument('--bcf_reference_panel', help='the reference panel bcf files [chr\tbcf]', required=True)
parser.add_argument('--genetic_map', help='the hg19 genetic map file', required=True)
parser.add_argument('--m3vcf_reference_panel', help='the reference panel m3vcf files [chr\tbcf]', required=True)

# DEFAULTS
parser.add_argument('--eagle', help='the path to eagle', required=False, default="eagle")
parser.add_argument('--minimac', help='the path to minimac3', required=False, default="minimac3")
parser.add_argument('--out_vcf_prefix', help='the prefix for the final vcf files', required=False, default="final")
parser.add_argument('--tabix', help='tabix executable location', required=False, default="tabix")
parser.add_argument('--eagle_threads', help='The number of threads for eagle', required=False, default=2, type=int)
parser.add_argument('--minimac_threads', help='he number of threads for minimac', required=False, default=2, type=int)
parser.add_argument('--delimiter', help='File list delimiter', required=False, default="\t")
args = parser.parse_args()

print

print "########################################################"
print "OPTIONS"
print "########################################################"

print

for attr, value in args.__dict__.iteritems():
        print "{0:>15}\t\t{1!s}".format( "--{0}".format(attr), str(value))

outdir = args.outdir
vcf_files_list = args.vcf_files
bcf_reference_panel_list = args.bcf_reference_panel
genetic_map = args.genetic_map
m3vcf_reference_panel_list = args.m3vcf_reference_panel
eagle = args.eagle
minimac = args.minimac
out_vcf_prefix = args.out_vcf_prefix
tabix = args.tabix
eagle_threads = args.eagle_threads
minimac_threads = args.minimac_threads
delimiter = args.delimiter

def read_list_file(file_path):
	map_of_files = dict()
	file_handle = open(file_path,'r')
	for line in file_handle:
		cols = line.rstrip().split(delimiter)
		map[cols[0]] = cols[1]
	file_handle.close()
	return map_of_files
	

vcf_files_map = read_list_file(vcf_files_list)
bcf_reference_map = read_list_file(bcf_reference_panel_list)
m3vcf_reference_map = read_list_file(m3vcf_reference_panel_list)


chroms = vcf_files_map.keys()

####
# Validate keys
####

for chrom in chroms:
	if chrom not in bcf_reference_map:
		raise ValueError("{chrom} not found in {bcfs}".format(chrom=chrom,bcfs=bcf_reference_panel_list))

	if chrom not in m3vcf_reference_map:
		raise ValueError("{chrom} not found in {m3s}".format(chrom=chrom,m3s=m3vcf_reference_panel_list))
		
eagle_vcfs = []
minimac_vcfs = []
final_vcfs = []
all_entries = []
for chrom in chroms:
	# EAGLE
	vcf = vcf_files_map[chrom]
	bcf_ref = bcf_reference_map[chrom]

	eagle_out_prefix = "{outdir}/eagle.phased.chr{chrom}".format(outdir=outdir,chrom=chrom)
	eagle_out_vcf = "{outdir}/eagle.phased.chr{chrom}.vcf.gz".format(outdir=outdir,chrom=chrom)
	eagle_vcfs.append( eagle_out_vcf  )
	
	eagle_cmd = "{eagle} --vcfTarget={vcf} --vcfRef={bcf_ref} --geneticMapFile={genetic_map} --outPrefix={eagle_out_prefix} --numThreads={eagle_threads} --vcfOutFormat=z".format(eagle=eagle,vcf=vcf,bcf_ref=bcf_ref,genetic_map=genetic_map,eagle_out_prefix=eagle_out_prefix,eagle_threads=eagle_threads)
	eagle_target = "{outdir}/eagle.phased.chr{chrom}.OK".format(outdir=outdir,chrom=chrom)
	eagle_entry = MakeEntry(eagle_target, [ eagle_cmd ], [ ], comment="Eagle phasing for chr{chrom}".format(chrom=chrom))
	all_entries.append(eagle_entry)
	
	# MINIMAC
	m3vcf = m3vcf_reference_map[chrom]
	minimac_out_prefix = "{outdir}/minimac.imputed.chr{chrom}".format(outdir=outdir,chrom=chrom)
	minimac_out_vcf = "{outdir}/minimac.imputed.chr{chrom}.dose.vcf.gz".format(outdir=outdir,chrom=chrom)
	minimac_cmd = "{minimac} --cpus {minimac_threads} --refHaps {m3vcf} --haps {eagle_out_vcf} --prefix {minimac_out_prefix}".format(minimac=minimac,minimac_threads=minimac_threads,m3vcf=m3vcf,eagle_out_vcf=eagle_out_vcf,minimac_out_prefix=minimac_out_prefix)
	minimac_target = "{outdir}/minimac.imputed.chr{chrom}.OK".format(outdir=outdir,chrom=chrom)
	minimac_entry = MakeEntry(minimac_target,  [ minimac_cmd ], [ eagle_entry ], comment="minimac phasing for chr{chrom}".format(chrom=chrom))
	all_entries.append(minimac_entry)
	
	#Rename final vcfs
	final_vcf = ""
	rename_cmd = "mv {vcf} {final_vcf}".format(vcf=minimac_out_vcf,final_vcf=final_vcf)
	tabix_cmd = "{tabix} -pvcf {final_vcf}".format(tabix=tabix,final_vcf=final_vcf)
	final_vcfs.append(final_vcf)
	rename_target = "{outdir}/final.chr{chrom}.OK".format(outdir=outdir,chrom=chrom)
	rename_entry = MakeEntry(rename_target, [rename_cmd,tabix_cmd], [minimac_entry], comment="Rename and finalize chr{chrom} vcf".format(chrom=chrom))
	
	