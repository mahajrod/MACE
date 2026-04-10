#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from RouToolPa.Parsers.VCF import CollectionVCF
from RouToolPa.Parsers.Sequence import CollectionSequence
from MACE.Routines import StatsVCF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with variants")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-r", "--reference_genome", action="store", dest="reference",
                    help="Fasta file with reference genome, required to get scaffold lengths."
                         "If absent lengths will be extracted from vcf metadata")
parser.add_argument("-w", "--window_size", action="store", dest="window_size", default=100000, type=int,
                    help="Size of the windows Default: 100000")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", default=None, type=int,
                    help="Step of the sliding windows. Default: window size, i.e windows are staking")
parser.add_argument("-e", "--per_sample", action="store_true", dest="per_sample", default=False,
                    help="Count variants for each sample independently. Default: count all samples together")

parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="parse",
                    help="Parsing mode for input sequence file. "
                         "Possible variants: 'index_db', 'index', 'parse'(default)")

args = parser.parse_args()

variants = CollectionVCF(in_file=args.input, parsing_mode='genotypes')

if args.reference:
    reference_length_df = CollectionSequence(in_file=args.reference, format="fasta", parsing_mode="parse", get_stats=True).seq_lengths
else:
    reference_length_df = None

StatsVCF.count_variants_in_windows(variants, args.window_size, args.window_step,
                                   reference_scaffold_lengths=reference_length_df,
                                   ignore_scaffolds_shorter_than_window=True, output_prefix=args.output_prefix,
                                   skip_empty_windows=False, expression=None, per_sample_output=args.per_sample)
