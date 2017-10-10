#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from BCBio import GFF
from MACE.Parsers.VCF import CollectionVCF, ReferenceGenome

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with variants.")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
"""
parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=300,
                    help="Dpi of figure")

parser.add_argument("-f", "--size_of_figure", action="store", dest="size_of_figure", type=lambda s: s.split(","),
                    default=(40, 40),
                    help="Size of figure in inches. X and Y values should be separated by comma. Default: 40,40")
"""
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=("svg", "png"),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg,png")

parser.add_argument("-l", "--suptitle", action="store", dest="suptitle", default="Variant density",
                    help="Suptitle of figure. Default: 'Variant density'")
"""
parser.add_argument("-g", "--draw_gaps", action="store_true", dest="draw_gaps",
                    help="Draw gaps, ignored if reference genome is not set. Default: False")
"""
parser.add_argument("-r", "--reference_genome", action="store", dest="reference",
                    help="Fasta file with reference genome, required to draw gaps and chromosomes")
"""
parser.add_argument("-m", "--masked_regions", action="store", dest="masked_regions",
                    help="Gff file with masked regions")
"""
parser.add_argument("-w", "--window_size", action="store", dest="window_size", default=100000, type=int,
                    help="Size of the windows Default: 100000")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", default=None, type=int,
                    help="Step of the sliding windows. Default: window size, i.e windows are staking")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="index_db",
                    help="Parsing mode for input sequence file. "
                         "Possible variants: 'index_db'(default), 'index', 'parse'")
args = parser.parse_args()

variants = CollectionVCF(from_file=True, in_file=args.input, parse_only_coordinates=True)
"""
if args.ref_genome:
    reference_genome = ReferenceGenome(args.reference)
    reference_genome.find_gaps()
else:
    reference_genome = None

if args.masked_regions:
    masked_regions = {}
    with open(args.masked_regions) as gff_fd:
        for record in GFF.parse(gff_fd):
            masked_regions[record.id] = record
else:
    masked_regions = None
"""

variants.draw_variant_window_densities(args.reference, args.output_prefix, args.window_size,
                                       args.window_size if args.window_step is None else args.window_step,
                                       masking=None, parsing_mode=args.parsing_mode, min_gap_length=10,
                                       masked_region_color="grey", gap_color="white",
                                       ignore_scaffolds_shorter_than_window=True,
                                       skip_empty_windows=False, scaffold_black_list=(),
                                       figure_extensions=args.output_formats)
