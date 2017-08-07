#!/usr/bin/python

import argparse
import re
import pysam
import os
import string
from pysam import VariantFile
import itertools

"""
annotate_gwas_results_from_vcf.py This script will take a vcf file and gwas results and add additional columns from a 
SnpEff annotated VCF
"""
