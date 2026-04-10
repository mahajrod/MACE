#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import VCFRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input vcf file with mutations")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    required=True,
                    help="Prefix of output files")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="one",
                    help="Operation mode. Allowed: 'one'(default) - variant will be treated as heterozygous if "
                         "there is at least one heterozygous sample, 'all' - all samples have to be heterozygous")

args = parser.parse_args()

VCFRoutines.extract_heterozygous_variants(args.input, args.output_prefix,
                                          mode=args.mode, verbose=True)
