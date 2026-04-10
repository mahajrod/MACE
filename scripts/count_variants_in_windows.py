#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from RouToolPa.Parsers.VCF import CollectionVCF
from RouToolPa.Parsers.Sequence import CollectionSequence
from MACE.Routines import StatsVCF
from MACE.Routines import Parsing

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with variants")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-n", "--length_file", action="store", dest="length_file",
                    help="File with sequence lengths. Allowed: '.fasta', '.fai', '.len' "
                         "File type is detected from the extension."
                         "If not set lengths will be extracted from vcf metadata")
parser.add_argument("-w", "--window_size", action="store", dest="window_size", default=100000, type=int,
                    help="Size of the windows Default: 100000")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", default=None, type=int,
                    help="Step of the sliding windows. Default: window size, i.e windows are staking")
parser.add_argument("-e", "--per_sample", action="store_true", dest="per_sample", default=False,
                    help="Count variants for each sample independently. Default: count all samples together")


args = parser.parse_args()

variants = CollectionVCF(in_file=args.input, parsing_mode='genotypes')

auxiliary_dict = Parsing.read_mace_auxiliary_input(len_file=args.length_file)

count_df = StatsVCF.count_variants_in_windows(variants, args.window_size, args.window_step,
                                              reference_scaffold_lengths=auxiliary_dict["len_df"],
                                              ignore_scaffolds_shorter_than_window=True, output_prefix=args.output_prefix,
                                              skip_empty_windows=False, expression=None, per_sample_output=args.per_sample)
feature_df, track_df = StatsVCF.convert_variant_count_to_feature_df(count_df, args.window_size,
                                                                    args.window_step if args.window_step is not None else args.window_size)

feature_df.to_csv("{}.counts.bed".format(args.output_prefix), sep="\t", header=True, index=True)
# TODO: refactor counting system later
track_columns = list(track_df.columns)
track_columns = track_columns[0:2] + [track_columns[-1]]
track_df[track_columns].to_csv("{}.counts.bedgraph".format(args.output_prefix), sep="\t", header=False, index=True)
