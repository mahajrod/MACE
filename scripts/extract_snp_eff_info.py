#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from MACE.Parsers.VCF import CollectionVCF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with mutations.")
parser.add_argument("-o", "--output_file", action="store", dest="output", required=True,
                    help="Output file with SNPeff info")
parser.add_argument("-e", "--snpeff_entry", action="store", dest="snpeff_entry", default="ANN",
                    help="SNPeff entry in vcf file. Default: ANN")

args = parser.parse_args()

mutations = CollectionVCF(in_file=args.input, from_file=True)
mutations.extract_snpeff_info(args.output, snpeff_entry=args.snpeff_entry)
