#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from MACE.Parsers.VCFpandas import CollectionVCF
from MACE.Parsers.Sequence import CollectionSequence

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with mutations.")
parser.add_argument("-o", "--output", action="store", dest="output",
                    required=True,
                    help="Output file with masking")
parser.add_argument("-s", "--sample_list", action="store", dest="sample_list",
                    default=None, type=lambda s: s.split(","),
                    help="Comma-separated list of samples to use for masking calculation. Default: all")
parser.add_argument("-u", "--min_samples", action="store", dest="min_samples",
                    default=1, type=int,
                    help="Minimum number of samples to mask variant. Default: 1")
parser.add_argument("-x", "--max_coverage", action="store", dest="max_coverage",
                    default=2.5, type=float,
                    help="Maximum coverage(relative to median) to retain variant. Default: 2.5")
parser.add_argument("-u", "--min_coverage", action="store", dest="min_coverage",
                    default=None, type=float,
                    help="Minimum coverage(relative to median) to retain variant. Default: not set")

args = parser.parse_args()

mutations = CollectionVCF(args.input, parsing_mode="pos_gt_dp")
mutations.calculate_masking(args.output, samples=args.sample_list, min_samples=args.min_samples,
                            max_coverage=args.max_coverage, min_coverage=args.min_coverage)
