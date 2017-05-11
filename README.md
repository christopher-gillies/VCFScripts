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
