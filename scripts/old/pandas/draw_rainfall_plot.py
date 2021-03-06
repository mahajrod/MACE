#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Parsers.VCF import CollectionVCF
from RouToolPa.Parsers.Sequence import CollectionSequence

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with mutations.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    required=True,
                    help="Prefix of output file with rainfall plot")
parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=300,
                    help="Dpi of figure")
parser.add_argument("-f", "--figsize", action="store", dest="figsize",
                    type=lambda s: map(int, s.split(",")),
                    default=(20, 20),
                    help="Size of figure in inches. X and Y values should be separated "
                         "by comma. Default: 20,20")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats",
                    type=lambda s: s.split(","),
                    default=["png"],
                    help="Comma-separated list of formats (supported by matlotlib) "
                         "of output figure.Default: png")
parser.add_argument("-l", "--suptitle", action="store", dest="suptitle",
                    help="Suptitle of figure. Default: 'Rainfall plot'")

parser.add_argument("-q", "--fontsize", action="store", dest="fontsize", type=int,
                    help="Size of font for suptitle and scaffold labels. Default: matplotlib defaults")

parser.add_argument("-r", "--reference_genome", action="store", dest="ref_genome",
                    help="Fasta file with reference genome, required to draw gaps")
parser.add_argument("-p", "--reference_parsing_mode", action="store",
                    dest="reference_parsing_mode", default="parse",
                    help="Parsing mode for fasta file with reference genome. Use 'generator' mode"
                         "if your machine have low memory. Allowed parse(default), generator")
parser.add_argument("-s", "--draw_masking", action="store_true", dest="draw_masking",
                    help="Draw masked regions. Default: false")
parser.add_argument("-m", "--masked_regions", action="store", dest="masked_regions",
                    help="Gff file with masked regions")
parser.add_argument("-n", "--min_masking_length", action="store", dest="min_masking_length", default=1,
                    type=int,
                    help="Minimum length of masked or gapped region to be shown. "
                         "Increase this parameter if you deal wit hlarge genome. Default: 1(show all)")
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
parser.add_argument("-w", "--dot_size", action="store", dest="dot_size", default=None, type=float,
                    help="Size of dots corresponding to the variants. Default: matplotlib default.")
args = parser.parse_args()

mutations = CollectionVCF(args.input, parsing_mode="only_coordinates")

if args.ref_genome:
    reference_genome = CollectionSequence(args.ref_genome,
                                          parsing_mode=args.reference_parsing_mode,
                                          masking_file=args.masked_regions,
                                          black_list=args.scaffold_black_list,
                                          white_list=args.scaffold_white_list)
    reference_genome.get_stats_and_features(count_gaps=True if args.draw_masking else False, sort="True", min_gap_length=1)
else:
    reference_genome = None

mutations.rainfall_plot(args.output_prefix, dpi=args.dpi, figsize=args.figsize,
                        facecolor="#D6D6D6",
                        ref_genome=reference_genome,
                        min_masking_length=args.min_masking_length,
                        suptitle=args.suptitle,
                        masking_color="#777777", logbase=args.logbase,
                        extension_list=args.output_formats,
                        scaffold_black_list=args.scaffold_black_list, scaffold_white_list=args.scaffold_white_list,
                        scaffold_ordered_list=args.scaffold_ordered_list, sort_scaffolds=args.sort_scaffolds,
                        color_expression=None,
                        default_point_color='blue',
                        dot_size=args.dot_size,
                        label_fontsize=args.fontsize,
                        draw_masking=args.draw_masking)
