#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from MACE.Parsers.SNPeff import CollectionSNPeff


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_snpeff", action="store", dest="input_snpeff", required=True,
                    help="Input snpeff tab file")
parser.add_argument("-o", "--output_snpeff", action="store", dest="output_snpeff", required=True,
                    help="Output snpeff file")

args = parser.parse_args()

snpeff_collection = CollectionSNPeff(input_file=args.input_snpeff, from_file=True, filetype='tab')
snpeff_collection.write(args.output_snpeff)


