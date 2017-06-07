#!/usr/bin/python

import argparse
import re
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
print "This script creates a makefile to concat vcfs with a different sample order"
print "########################################################"
parser = argparse.ArgumentParser(description='This program takes a plink genotype file and formats it for the Imputation Server')
parser.add_argument('--outdir', help='Output file prefix', required=True)
parser.add_argument('--vcf_list', help='List of vcf files to concat [chr tab vcf]', required=True)
parser.add_argument('--sample_list', help='The samples to use', required=True)
parser.add_argument('--maf', help='minimum minor allele frequency', type=float, required=False, default=0.05)
parser.add_argument('--threads', help='threads', required=False, type=int, default=2)
args = parser.parse_args()

print

print "########################################################"
print "OPTIONS"
print "########################################################"

print

for attr, value in args.__dict__.iteritems():
        print "{0:>15}\t\t{1!s}".format( "--{0}".format(attr), str(value))

vcf_list = args.vcf_list
sample_list = args.sample_list
maf = args.maf
outdir = args.outdir
min_af = maf
max_af = 1 - maf
threads = args.threads


print
print "CREATING MAKEFILE"
print 

make_entries = []

vcfs_to_clear = []
reorder_entries = []
with open(vcf_list,"r") as vcf_list_h:
	for line in vcf_list_h:
		line = line.rstrip()
		chrom,vcf = line.split("\t")
		match = re.search("(chr)?(.+)",chrom)
		chrom_num = match.group(2)
		out_vcf = "{0}/chr{1}.bcf".format(outdir,chrom_num)
		vcfs_to_clear.append(out_vcf)
	
		cmd = "bcftools view {0} --samples-file {1} --output-type b --output-file {2} --min-af {3} --max-af {4}".format(vcf,sample_list,out_vcf,min_af,max_af)
		reorder_entry = MakeEntry(  "{0}/reorder.chr{1}.OK".format(outdir,chrom_num), [ cmd ], [ ], comment='Reorder chr{0}'.format(chrom_num))
		make_entries.append(reorder_entry)
		reorder_entries.append(reorder_entry)

bcf_list = "{0}/bcfs".format(outdir)

with open(bcf_list,"w") as bcf_h:
	for vcf in vcfs_to_clear:
		bcf_h.write(vcf)
		bcf_h.write("\n")

concat_cmd = "bcftools concat -f {0} --output-type b --output {1}/merged.bcf --threads {2}".format(bcf_list,outdir,threads)
index_cmd = "bcftools index {0}/merged.bcf".format(outdir)

concat_entry = MakeEntry(  "{0}/concat.OK".format(outdir), [ concat_cmd, index_cmd ], reorder_entries , comment='Concat command')
make_entries.append(concat_entry)

clean_cmds = []
for vcf in vcfs_to_clear:
	cmd = "rm {0}".format(vcf)
	clean_cmds.append(cmd)

clean_entry = MakeEntry("clean", clean_cmds, [], comment ="Clean command")


makefile = MakeFile("{0}/all.OK".format(outdir), make_entries)
makefile.set_clean_entry(clean_entry)
makefile_loc = "{0}/concat.make".format(outdir)
out_makefile = open(makefile_loc, "w")
out_makefile.write(str(makefile))
out_makefile.write("\n")
out_makefile.close()

print "MAKEFILE CREATED"
print "Please run:"
print "cd " + outdir + ";" + "make -f concat.make"
