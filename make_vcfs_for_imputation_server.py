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
print "Genotype File formatting pipeline for imputation server"
print "########################################################"
parser = argparse.ArgumentParser(description='This program takes a plink genotype file and formats it for the Imputation Server')
parser.add_argument('--bfile', help='plink formatted binary genotype file prefix (same as plink)', required=True)
parser.add_argument('--outdir', help='the output directory ', required=True)
parser.add_argument('--strandfile', help='Download from http://www.well.ox.ac.uk/~wrayner/strand/', required=True)
parser.add_argument('--checkVCF', help='checkVCF https://github.com/zhanxw/checkVCF', required=True)
parser.add_argument('--refseq', help='Reference fasta file with checkVCF.py', required=True)
parser.add_argument('--hrcsites', help='Download from ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz', required=True)
parser.add_argument('--hrccheck', help='Downlaod from http://www.well.ox.ac.uk/~wrayner/tools/#Checking', required=True)

# DEFAULTS
parser.add_argument('--out_vcf_prefix', help='the prefix for the final vcf files', required=False, default="final")
parser.add_argument('--sort_buffer', help='Buffer size for sort program', required=False, default="4G")
parser.add_argument('--plink', help='plink executable location', required=False, default="plink")
parser.add_argument('--plink_threads', help='Number of threads for plink', required=False, default=1, type=int)
parser.add_argument('--tabix', help='tabix executable location', required=False, default="tabix")
parser.add_argument('--flip_ref_alt', help='flip_ref_alt.py executable location', required=False, default="{0}/{1}".format(calling_dir,"flip_ref_alt.py"))
parser.add_argument('--mind', help='plink --mind parameter (default 0.05)', type=float, default=0.05)
parser.add_argument('--geno', help='plink --geno parameter (default 0.05)', type=float, default=0.05)
args = parser.parse_args()

print

print "########################################################"
print "OPTIONS"
print "########################################################"

print

for attr, value in args.__dict__.iteritems():
        print "{0:>15}\t\t{1!s}".format( "--{0}".format(attr), str(value))

##############
# Store inputs
##############
bfile = args.bfile
outdir = args.outdir

if outdir.endswith("/"):
	outdir = outdir.rstrip("/")
	
outdir_temp = "{0}/tmp".format(outdir)
strandfile = args.strandfile
checkVCF = args.checkVCF
refseq = args.refseq
hrcsites = args.hrcsites
hrccheck = args.hrccheck

out_vcf_prefix = args.out_vcf_prefix
sort_buffer = args.sort_buffer
plink = args.plink
plink_threads = args.plink_threads
#vcfsort = args.vcfsort
tabix = args.tabix
flip_ref_alt = args.flip_ref_alt
mind = args.mind
geno = args.geno

#mk_test = MakeEntry("11111", ["print a","print b"], [])	
#mk_test_a = MakeEntry("2222", ["print c","print d"], [mk_test])
#mk_file = MakeFile("all.ok", [ mk_test,mk_test_a ])
#print mk_file


print
print "CREATING MAKEFILE"
print 

make_entries = []
############
# Create Commands
############

########
# Update positions, extact snps in strand file, filter geno/inds, flip strands, make frequency file
########

bfile_base = os.path.basename(bfile)
make_outdir_tmp_cmd = "mkdir -p {0}".format(outdir_temp)

snp_id_file = "{0}/snp.ids".format(outdir_temp)
snp_id_command = "cat {0} | cut -f1 > {1}".format(strandfile,snp_id_file)

minus_strand_file = "{0}/minus.snp.ids".format(outdir_temp)
minus_strand_cmd = "cat {strandfile} | perl -lane 'print $$F[0] if $$F[4] eq \"-\"' > {minus_strand_file}".format(strandfile=strandfile, minus_strand_file=minus_strand_file)

bfile_plink_geno_prefix = "{0}/extract.geno".format(outdir_temp)
extract_snps_geno_filter_cmd = "{plink} --bfile {bfile} --geno {geno} --extract {snp_id_file} --make-bed --out {out} --threads {plink_threads}".format(plink=plink, bfile=bfile, geno=geno, snp_id_file=snp_id_file, out=bfile_plink_geno_prefix,plink_threads=plink_threads)

bfile_plink_geno_pos_ind_prefix = "{0}/extract.geno.pos.ind".format(outdir_temp)
update_pos_ind_cmd = "{plink} --bfile {bfile} --mind {mind} --update-map {strandfile} 3 1 --make-bed --out {out} --threads {plink_threads}".format(plink=plink, bfile=bfile_plink_geno_prefix, strandfile=strandfile, mind=mind, out=bfile_plink_geno_pos_ind_prefix, plink_threads=plink_threads)

bfile_plink_geno_pos_ind_flip_prefix = "{0}/extract.geno.pos.ind.flip".format(outdir_temp)
flip_cmd = "{plink} --bfile {bfile} --flip {minus_strand_file} --make-bed --freq --out {out} --threads {plink_threads}".format(plink=plink, bfile=bfile_plink_geno_prefix, strandfile=strandfile, minus_strand_file=minus_strand_file, out=bfile_plink_geno_pos_ind_flip_prefix,plink_threads=plink_threads)


pre_hrc_check_entry = MakeEntry(  "{0}/pre_hrc_check.OK".format(outdir_temp), [make_outdir_tmp_cmd, snp_id_command, minus_strand_cmd, extract_snps_geno_filter_cmd, update_pos_ind_cmd, flip_cmd], [], comment = "PRE-HRC CHECKING FORMATING" )

flip_freq_file = "{0}/extract.geno.pos.ind.flip.frq".format(outdir_temp)
flip_bim_file = "{0}/extract.geno.pos.ind.flip.bim".format(outdir_temp)

make_entries.append(pre_hrc_check_entry)

#######
# HRC-Check
#######

hrc_cmd = "cd {outdir_temp}; {hrccheck} --bim {flip_bim_file} --freq {flip_freq_file} --ref {hrcsites} --hrc -v".format(outdir_temp=outdir_temp,hrccheck=hrccheck, flip_bim_file=flip_bim_file, flip_freq_file=flip_freq_file, hrcsites=hrcsites)
hrc_check_entry = MakeEntry(  "{0}/hrc_check.OK".format(outdir_temp), [hrc_cmd], [pre_hrc_check_entry], comment = "HRC CHECK " )

make_entries.append(hrc_check_entry)

#######
# HRC FORMAT
#######

#append plink threads
hrc_format_script = "cat {outdir_temp}/Run-plink.sh | perl -lane 'print $$_.\" --threads {plink_threads} --allow-extra-chr \" if $$_ =~ /^plink/'  > {outdir_temp}/Run-plink-updated.sh".format(outdir_temp=outdir_temp,plink_threads=plink_threads)
hrc_run_cmd = "cd {outdir_temp}; sh {outdir_temp}/Run-plink-updated.sh".format(outdir_temp=outdir_temp)
hrc_format_entry = MakeEntry(  "{0}/hrc_format.OK".format(outdir_temp), [hrc_format_script , hrc_run_cmd], [hrc_check_entry], comment = "HRC FORMAT" )
make_entries.append(hrc_format_entry)

#######
# Convert Plink Files to VCF sort and tabix
#######

# {outdir_temp}/extract.geno.pos.ind.flip-updated-chr7
hrc_bfile_prefix = []
hrc_vcf_file_prefix = []
hrc_vcf_file_out = []
hrc_vcf_file_out_sorted = []
hrc_vcf_file_out_sorted_prefix = []
vcf_convert_targets = []
vcf_sort_targets = []

# Chr start at 1
for i in range(1,24):
	bfile_path = "{0}/extract.geno.pos.ind.flip-updated-chr{1}".format(outdir_temp,i)
	vcf_path = "{0}/extract.geno.pos.ind.flip.updated.chr{1}.vcf.gz".format(outdir_temp,i)
	vcf_sorted_path = "{0}/extract.geno.pos.ind.flip.updated.chr{1}.sorted.vcf.gz".format(outdir_temp,i)
	vcf_sorted_prefix = "{0}/extract.geno.pos.ind.flip.updated.chr{1}.sorted".format(outdir_temp,i)
	vcf_convert = "{0}/vcf.convert.chr{1}.OK".format(outdir_temp,i)
	vcf_sort_tar = "{0}/vcf.sort.chr{1}.OK".format(outdir_temp,i)
	vcf_prefix = "{0}/extract.geno.pos.ind.flip.updated.chr{1}".format(outdir_temp,i)
	#hrc_bfile_prefix.append( os.path.splitext(bfile_path)[0]  )
	hrc_bfile_prefix.append( bfile_path )
	hrc_vcf_file_prefix.append(vcf_prefix)
	hrc_vcf_file_out.append(vcf_path)
	hrc_vcf_file_out_sorted.append(vcf_sorted_path)
	hrc_vcf_file_out_sorted_prefix.append(vcf_sorted_prefix)
	vcf_convert_targets.append(vcf_convert)
	vcf_sort_targets.append(vcf_sort_tar)

####
# Create Sort and Convert commands
####
convert_entries = []
sort_entries = []
# Array Start at 0	
for i in range(0,23):
	plink_to_vcf_cmd = "{plink} --bfile {bfile} --recode vcf bgz --out {vcf_prefix} --threads {plink_threads}".format(plink=plink,bfile=hrc_bfile_prefix[i],vcf_prefix=hrc_vcf_file_prefix[i], plink_threads=plink_threads)
	#vcf_sort_cmd = "{vcfsort} {vcf} -t {tmp} | bgzip -c > {sorted_vcf}".format(vcfsort=vcfsort,vcf=hrc_vcf_file_out[i],sorted_vcf=hrc_vcf_file_out_sorted[i],tmp=outdir_temp)
	
	vcf_header_cmd = "zcat {vcf} | head -10000 | grep '^#' > {sorted_vcf}.header".format(vcf=hrc_vcf_file_out[i],sorted_vcf=hrc_vcf_file_out_sorted[i])
	
	#vcf_sort_cmd = "( cat {sorted_vcf}.header; (zcat {vcf} | grep -v '^#' | sort -t$$'\\t' -T {tmp} -S {sort_buffer} -k1,1 -k2,2n) ) | bgzip -c > {sorted_vcf}".format(sort_buffer=sort_buffer,vcf=hrc_vcf_file_out[i],sorted_vcf=hrc_vcf_file_out_sorted[i],tmp=outdir_temp)
	vcf_sort_cmd = "( cat {sorted_vcf}.header; (zcat {vcf} | grep -v '^#' | sort -t$$'\\t' -S {sort_buffer} -k1,1 -k2,2n) ) | bgzip -c > {sorted_vcf}".format(sort_buffer=sort_buffer,vcf=hrc_vcf_file_out[i],sorted_vcf=hrc_vcf_file_out_sorted[i],tmp=outdir_temp)
	
	tabix_cmd = "{tabix} -pvcf {vcf}".format(tabix=tabix,vcf=hrc_vcf_file_out_sorted[i])
	
	convert_to_vcf_entry = MakeEntry(vcf_convert_targets[i], [ plink_to_vcf_cmd ], [ hrc_format_entry ], comment="VCF conversion from plink for chr{0}".format(i + 1))
	sort_entry = MakeEntry(vcf_sort_targets[i], [ vcf_header_cmd, vcf_sort_cmd, tabix_cmd ], [ convert_to_vcf_entry ], comment="Sort VCF from plink for chr{0}".format(i + 1))
	
	make_entries.append(convert_to_vcf_entry)
	make_entries.append(sort_entry)
	
	convert_entries.append(convert_to_vcf_entry)
	sort_entries.append(sort_entry)
	
#######
# run checkVCF.py
#######
check_vcf_out = []
check_vcf_out_ref_flips = []
check_vcf_target = []
check_vcf_entries = []

for i in range(0,23):
	vcf = hrc_vcf_file_out_sorted[i]
	target = "{0}/checkVCF.chr{1}.OK".format(outdir_temp,i + 1)
	check_vcf_out.append(check_vcf_out)
	check_vcf_target.append(target)
	prefix = hrc_vcf_file_out_sorted_prefix[i]
	check_vef_ref_flips = "{0}.check.ref".format(prefix)
	check_vcf_out_ref_flips.append(check_vef_ref_flips)
	check_vcf_cmd = "python {checkVCF} -r {refseq} -o {prefix} {vcf}".format(checkVCF =checkVCF, refseq = refseq, prefix = prefix, vcf = vcf)
	check_entry = MakeEntry(target, [ check_vcf_cmd ], [ convert_entries[i], sort_entries[i] ], comment = "checkVCF for chr{0}".format(i + 1))
	make_entries.append(check_entry)
	check_vcf_entries.append(check_entry)

#######
# run flip remaining variants
#######
flip_ref_alt_vcf_out = []
flip_ref_alt_entries = []
for i in range(0,23):
	in_vcf = hrc_vcf_file_out_sorted[i]
	target = "{0}/flip_ref_alt.chr{1}.OK".format(outdir_temp,i + 1)
	out_vcf = "{0}/extract.geno.pos.ind.flip.updated.chr{1}.sorted.fix_ref_alts.vcf.gz".format(outdir,i + 1)
	flip_ref_alt_vcf_out.append(out_vcf)
	flip_ref_alt_cmd = "{flip_ref_alt} --vcf {in_vcf} --ref_flips {flips} | bgzip -c > {out_vcf}".format(flip_ref_alt=flip_ref_alt,in_vcf=in_vcf,out_vcf=out_vcf,flips=check_vcf_out_ref_flips[i])	
	tabix_cmd = "{tabix} -pvcf {out_vcf}".format(tabix=tabix, out_vcf=out_vcf)
	flip_ref_alt_entry = MakeEntry(target, [ flip_ref_alt_cmd, tabix_cmd], [ check_vcf_entries[i] ], comment = "flip ref and alt for chr{0}".format(i + 1))
	flip_ref_alt_entries.append(flip_ref_alt_entry)
	make_entries.append(flip_ref_alt_entry)


final_vcfs = []
for i in range(0,23):
	target = "{0}/final.vcf.chr{1}.OK".format(outdir_temp,i + 1)
	in_vcf = flip_ref_alt_vcf_out[i]
	out_vcf = "{outdir}/{out_vcf_prefix}.chr{chr}.vcf.gz".format( outdir=outdir, out_vcf_prefix=out_vcf_prefix, chr=i + 1)
	final_vcfs.append(out_vcf)
	final_cmd = "mv {in_vcf} {out_vcf}; mv {in_vcf}.tbi {out_vcf}.tbi".format(in_vcf=in_vcf,out_vcf=out_vcf)
	final_entry = MakeEntry(target, [final_cmd], [ flip_ref_alt_entries[i] ], comment="Rename final vcfs")
	make_entries.append(final_entry)
	
######
######

#########
# Clean-up
#########

clean_cmd = "rm -r {outdir_temp}".format(outdir_temp=outdir_temp)
clean_entry = MakeEntry("clean", [ clean_cmd ], [], comment ="Clean command")

###
# write out final vcf list
###
out_vcf_file_list = open( "{0}/vcf_list.tsv".format(outdir), "w")
for i in range(0,23):
	out_vcf_file_list.write("chr{chr}\t{vcf}\n".format(chr=i+1,vcf=final_vcfs[i]))
out_vcf_file_list.close()
	
####
# Create Makefile
####

makefile = MakeFile("{0}/all.OK".format(outdir), make_entries)
makefile.set_clean_entry(clean_entry)
makefile_loc = "{0}/Makefile".format(outdir)
out_makefile = open(makefile_loc, "w")
out_makefile.write(str(makefile))
out_makefile.write("\n")
out_makefile.close()

print "MAKEFILE CREATED"
print "Please run:"
print "cd " + outdir + ";" + "make"
