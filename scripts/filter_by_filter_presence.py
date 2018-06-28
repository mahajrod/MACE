#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from MACE.Parsers.VCF import CollectionVCF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input vcf file with variants")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

variants = CollectionVCF(in_file=args.input, from_file=True)
without_filters, with_filters = variants.filter_by_filter_presence()

without_filters.write("%s.without_filters.vcf" % args.output_prefix)
with_filters.write("%s.with_filters.vcf" % args.output_prefix)
