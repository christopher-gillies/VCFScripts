# VCFScripts
A set of simple VCF parsing scripts


# flip_ref_alt.py
run checkVCF on an input vcf to find out what needs to be swapped
https://github.com/zhanxw/checkVCF
the ref_flips input is .check.ref file 
```
export vcf=/path/to/vcf
export ref_flips=/path/to
/path/to/flip_ref_alt.py --vcf /path/to/vcf --ref_flips /path/to/ref_flips/ | bgzip -c > /path/to/out.vcf.gz
tabix /path/to/out.vcf.gz
```
# illumina_to_plink_lgen.py
## Convert from Illumina report to plink LGEN and BED format 

```
export format=/path/illumina_to_plink_lgen.py
export manifest=/path/manifest.csv
export sample_file=/path/sample_illumina_reports.tsv
$format --out_prefix samples.out --sample_file $sample_file --manifest_csv $manifest --skip_indels
plink --lfile samples.out  --make-bed --out samples.out  
```
