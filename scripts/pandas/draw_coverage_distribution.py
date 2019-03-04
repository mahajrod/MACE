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
parser.add_argument("-s", "--subplot_size", action="store", dest="subplot_size",
                    type=int, default=3,
                    help="Size of figure per sample subplot in inches. Default: 3")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats",
                    type=lambda s: s.split(","),
                    default=["png"],
                    help="Comma-separated list of formats (supported by matlotlib) "
                         "of output figure.Default: png")
parser.add_argument("-w", "--bin_width", action="store", dest="bin_width",
                    default=5, type=int,
                    help="Bin width for distribution histogram. Default: 5")
parser.add_argument("--verbose", action="store_true", dest="verbose",
                    help="Verbose stdout")
"""
parser.add_argument("-a", "--scaffold_white_list", action="store", dest="scaffold_white_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of the only scaffolds to draw. Default: all")
parser.add_argument("-b", "--scaffold_black_list", action="store", dest="scaffold_black_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of scaffolds to skip at drawing. Default: not set")
"""
args = parser.parse_args()

mutations = CollectionVCF(args.input, parsing_mode="pos_gt_dp")
mutations.get_coverage_distribution(args.output_prefix,
                                    bin_width=args.bin_width, dpi=args.dpi, subplot_size=3,
                                    extension_list=args.output_formats, verbose=args.verbose)
