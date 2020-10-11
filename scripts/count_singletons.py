#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Parsers.VCF import CollectionVCF
from MACE.Routines import StatsVCF, Visualization

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
                         "by comma. Default: 5,5")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats",
                    type=lambda s: s.split(","),
                    default=["png"],
                    help="Comma-separated list of formats (supported by matlotlib) "
                         "of output figure.Default: png")
parser.add_argument("-l", "--title", action="store", dest="title",
                    default=None,
                    help="Title of figure. Default: not set")
parser.add_argument("-m", "--parsing_mode", action="store", dest="parsing_mode",
                    default="genotypes",
                    help="Parsing mode. Allowed: genotypes(default), 'coordinates_and_genotypes', 'complete'")
"""
parser.add_argument("-a", "--scaffold_white_list", action="store", dest="scaffold_white_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of the only scaffolds to draw. Default: all")
parser.add_argument("-b", "--scaffold_black_list", action="store", dest="scaffold_black_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of scaffolds to skip at drawing. Default: not set")
"""
args = parser.parse_args()

mutations = CollectionVCF(args.input, parsing_mode="genotypes")
StatsVCF.count_singletons(collection_vcf=mutations, output_prefix=args.output_prefix)


"""
Visualization.zygoty_bar_plot(StatsVCF.count_zygoty(mutations, outfile="%s.counts" % args.output_prefix),
                              args.output_prefix, extension_list=args.output_formats,
                              figsize=args.figsize,
                              dpi=args.dpi,
                              title=args.title)
"""