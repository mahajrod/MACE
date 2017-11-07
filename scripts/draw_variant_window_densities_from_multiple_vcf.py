#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from collections import OrderedDict
from BCBio import GFF
from MACE.Parsers.VCF import CollectionVCF, ReferenceGenome
from MACE.Routines import DrawingRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of vcf file with variants.")
parser.add_argument("-n", "--sample_names", action="store", dest="sample_names", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of sample names with variants. They must be in same order as vcf files.")
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

args = parser.parse_args()

count_dict = OrderedDict()

reference = ReferenceGenome(args.reference,
                            masked_regions=None,
                            index_file="refgen.idx",
                            filetype="fasta",
                            mode=args.parsing_mode,
                            black_list=[])

reference.find_gaps(min_gap_length=10)

for sample_name, vcf_file in zip(args.sample_names, args.input):
    count_dict[sample_name] = CollectionVCF(from_file=True, in_file=vcf_file, parse_only_coordinates=True).count_variants_in_windows(args.window_size,
                                                                                                                                     args.window_size if args.window_step is None else args.window_step,
                                                                                                                                     reference.region_length,
                                                                                                                                     ignore_scaffolds_shorter_than_window=True,
                                                                                                                                     output_prefix="%s.%s" % (args.output_prefix, sample_name),
                                                                                                                                     skip_empty_windows=False)


DrawingRoutines.draw_variant_window_densities(count_dict, reference.region_length, args.window_size,
                                              args.window_size if args.window_step is None else args.window_step,
                                              args.output_prefix,
                                              record_style=None, ext_list=args.output_formats,
                                              label_fontsize=13, left_offset=0.2, figure_width=8,
                                              scaffold_synonym_dict=None,
                                              id_replacement_mode="partial", suptitle=None, density_multiplicator=1000,
                                              scaffold_black_list=args.scaffold_black_list,
                                              sort_scaffolds=args.sort_scaffolds,
                                              scaffold_ordered_list=args.scaffold_ordered_list,
                                              scaffold_white_list=args.scaffold_white_list,
                                              add_sample_name_to_labels=True,
                                              gap_color="grey",
                                              colormap_tuple_list=((0.0, "#333a97"), (0.1, "#3d3795"),
                                                                   (0.5, "#5d3393"), (0.75, "#813193"),
                                                                   (1.0, "#9d2d7f"), (1.25, "#b82861"),
                                                                   (1.5, "#d33845"), (2.0, "#ea2e2e"),
                                                                   (2.5, "#f5ae27")))




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

"""
variants.draw_variant_window_densities(args.reference, args.output_prefix, args.window_size,
                                       args.window_size if args.window_step is None else args.window_step,
                                       masking=None, parsing_mode=args.parsing_mode, min_gap_length=10,
                                       masked_region_color="grey", gap_color="white",
                                       ignore_scaffolds_shorter_than_window=True,
                                       skip_empty_windows=False, scaffold_black_list=args.scaffold_black_list,
                                       sort_scaffolds=args.sort_scaffolds,
                                       scaffold_ordered_list=args.scaffold_ordered_list,
                                       scaffold_white_list=args.scaffold_white_list,
                                       figure_extensions=args.output_formats,
                                       add_sample_name_to_labels=False,
                                       sample_label="SampleZZZ",)
"""
