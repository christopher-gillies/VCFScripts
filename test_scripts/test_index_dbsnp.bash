

export dbsnp=/kidneyomics/references/dbsnp_138.hg19.relabel.vcf.gz
export extract_haplotypes=/nfs/turbo/wgsa/cgillies/VCFScripts/extract_haplotypes.py 
export dbsnp_out=/kidneyomics/references/dbsnp_138.hg19.relabel.pickle

$extract_haplotypes --vcf_file $dbsnp --index_dbsnp --out_file $dbsnp_out

export vcf_file=/nfs/turbo/wgsa/cgillies/dbgap/CKID/processed_data/phs000524/phs000524_cg1_strand_hg19_forward_updated_chr6_sorted_ref_flip_phased_impute_1kg.dose.vcf.gz 
$extract_haplotypes --vcf_file $vcf_file --dbsnp $dbsnp_out --out_file $dbsnp_out --marker_ids rs113529773

