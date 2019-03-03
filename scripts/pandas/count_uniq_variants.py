#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from MACE.Parsers.VCFpandas import CollectionVCF
from MACE.Parsers.Sequence import CollectionSequence

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with mutations.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    required=True,
                    help="Prefix of output files")
parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=200,
                    help="Dpi of figure")
parser.add_argument("-f", "--figsize", action="store", dest="figsize",
                    type=lambda s: map(int, s.split(",")),
                    default=(5, 5),
                    help="Size of figure in inches. X and Y values should be separated "
                         "by comma. Default: 40,40")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats",
                    type=lambda s: s.split(","),
                    default=["png"],
                    help="Comma-separated list of formats (supported by matlotlib) "
                         "of output figure.Default: png")
parser.add_argument("-l", "--title", action="store", dest="title",
                    default="Unique variants",
                    help="Title of figure. Default: Unique variants")

"""
parser.add_argument("-a", "--scaffold_white_list", action="store", dest="scaffold_white_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of the only scaffolds to draw. Default: all")
parser.add_argument("-b", "--scaffold_black_list", action="store", dest="scaffold_black_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of scaffolds to skip at drawing. Default: not set")
"""
args = parser.parse_args()

mutations = CollectionVCF(args.input, parsing_mode="coordinates_and_genotypes")
mutations.count_uniq_variants(args.output_prefix,
                              extension_list=args.output_formats,
                              figsize=args.figsize,
                              dpi=args.dpi,
                              title=args.title)
