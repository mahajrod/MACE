#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from BCBio import GFF
from MACE.Parsers.VCF import CollectionVCF, ReferenceGenome

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input fasta file with reference.")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-m", "--masked_regions", action="store", dest="masked_regions", type=lambda s: s.split(","),
                    help="Comma-separated list of Gff file with masked regions")

parser.add_argument("-w", "--window_size", action="store", dest="window_size", default=100000, type=int,
                    help="Size of the windows Default: 100000")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", default=None, type=int,
                    help="Step of the sliding windows. Default: window size, i.e windows are staking")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="parse",
                    help="Parsing mode for input sequence file. "
                         "Possible variants: 'index_db', 'index', 'parse'(default)")


args = parser.parse_args()

reference = ReferenceGenome(args.reference,
                            masked_regions=None,
                            index_file="refgen.idx",
                            filetype="fasta",
                            mode=args.parsing_mode,
                            black_list=[],
                            masking_gff_list = args.masked_regions)

reference.count_gaped_and_masked_positions_in_windows(args.window_size, args.window_step,
                                                      ignore_scaffolds_shorter_than_window=True,
                                                      output_prefix=args.output_prefix,
                                                      min_gap_len=1)




