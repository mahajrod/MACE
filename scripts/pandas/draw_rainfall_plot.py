#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from BCBio import GFF
from MACE.Parsers.VCFpandas import CollectionVCF, ReferenceGenome

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with mutations.")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix",
                    required=True,
                    help="Prefix of output file with rainfall plot")
parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=300,
                    help="Dpi of figure")
parser.add_argument("-f", "--size_of_figure", action="store", dest="size_of_figure",
                    type=lambda s: map(int, s.split(",")),
                    default=(20, 20),
                    help="Size of figure in inches. X and Y values should be separated "
                         "by comma. Default: 40,40")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats",
                    type=lambda s: s.split(","),
                    default=["png"],
                    help="Comma-separated list of formats (supported by matlotlib) "
                         "of output figure.Default: png")
parser.add_argument("-l", "--suptitle", action="store", dest="suptitle",
                    help="Suptitle of figure. Default: 'Rainfall plot'")
parser.add_argument("-g", "--draw_gaps", action="store_true", dest="draw_gaps",
                    help="Draw gaps, ignored if reference genome is not set. Default: False")
parser.add_argument("-r", "--reference_genome", action="store", dest="ref_genome",
                    help="Fasta file with reference genome, required to draw gaps")
parser.add_argument("-m", "--masked_regions", action="store", dest="masked_regions",
                    help="Gff file with masked regions")
parser.add_argument("-u", "--logbase", action="store", dest="logbase", default=2, type=int,
                    help="Logbase of y axis")
parser.add_argument("-a", "--scaffold_white_list", action="store", dest="scaffold_white_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of the only scaffolds to draw. Default: all")
parser.add_argument("-b", "--scaffold_black_list", action="store", dest="scaffold_black_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of scaffolds to skip at drawing. Default: not set")
parser.add_argument("-y", "--sort_scaffolds", action="store_true", dest="sort_scaffolds", default=False,
                    help="Order  scaffolds according to their names. Default: False")
parser.add_argument("-z", "--scaffold_ordered_list", action="store", dest="scaffold_ordered_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of scaffolds to draw first and exactly in same order. "
                         "Scaffolds absent in this list are drawn last and in order according to vcf file . "
                         "Default: not set")
parser.add_argument("-w", "--dot_size", action="store", dest="dot_size", default=None, type=int,
                    help="Size of dots corresponding to the variants. Default: matplotlib default.")
args = parser.parse_args()

mutations = CollectionVCF(args.input, dont_parse_info_and_data=True)

if args.ref_genome:
    reference_genome = ReferenceGenome(args.ref_genome)
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

mutations.rainfall_plot(args.output_prefix, single_fig=True, dpi=args.dpi, figsize=args.size_of_figure,
                        facecolor="#D6D6D6",
                        ref_genome=reference_genome, masked_scaffolds=masked_regions, min_gap_length=10,
                        draw_gaps=args.draw_gaps, suptitle=args.suptitle,
                        gaps_color="#777777", masked_scaffolds_color="#aaaaaa", logbase=args.logbase,
                        extension_list=args.output_formats,
                        scaffold_black_list=args.scaffold_black_list, scaffold_white_list=args.scaffold_white_list,
                        scaffold_ordered_list=args.scaffold_ordered_list, sort_scaffolds=args.sort_scaffolds,
                        color_expression=None,
                        default_point_color='black',
                        dot_size=args.dot_size,
                        label_fontsize=None)
