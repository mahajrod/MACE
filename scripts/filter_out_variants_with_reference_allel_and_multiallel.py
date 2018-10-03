#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from MACE.Parsers.VCF import CollectionVCF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input vcf file with variants")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-m", "--max_allels", action="store", dest="max_allels", type=int,
                    help="Max allels per variant. Default: not set")
args = parser.parse_args()

variants = CollectionVCF(in_file=args.input, from_file=True)
without_filters, with_filters = variants.filter_variants_with_reference_allel_and_multiallelic(sample_index=None,
                                                                                               max_allels=None)

without_filters.write("%s.no_reference_allel%s.vcf" % (args.output_prefix,
                                                       ".max_allels_%i" % args.max_allels if args.max_allels else ""))
with_filters.write("%s.with_reference_allel%s.vcf" % (args.output_prefix,
                                                      ".more_then_allels_%i" % args.max_allels if args.max_allels else ""))
