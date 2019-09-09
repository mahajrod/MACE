#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Parsers.VCF import CollectionVCF


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input vcf file with mutations")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    required=True,
                    help="Prefix of output files")

args = parser.parse_args()

mutations = CollectionVCF(args.input, parsing_mode="all_no_parsing")

mutations.write(args.output_prefix, format='vcf', samples=None, split_samples=True)
