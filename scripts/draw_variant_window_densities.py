#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

import pandas as pd
from BCBio import GFF
from RouToolPa.Collections.General import SynDict, IdList
from RouToolPa.Parsers.VCF import CollectionVCF
from MACE.Routines import Visualization, StatsVCF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input vcf file with variants.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
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

parser.add_argument("-l", "--title", action="store", dest="title", default="Variant density",
                    help="Suptitle of figure. Default: 'Variant density'")
"""
parser.add_argument("-g", "--draw_gaps", action="store_true", dest="draw_gaps",
                    help="Draw gaps, ignored if reference genome is not set. Default: False")

parser.add_argument("-r", "--reference_genome", action="store", dest="reference",
                    help="Fasta file with reference genome, required to draw gaps and chromosomes")

parser.add_argument("-m", "--masked_regions", action="store", dest="masked_regions",
                    help="Gff file with masked regions")
"""
parser.add_argument("-w", "--window_size", action="store", dest="window_size", default=100000, type=int,
                    help="Size of the windows Default: 100000")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", default=None, type=int,
                    help="Step of the sliding windows. Default: window size, i.e windows are staking")

parser.add_argument("-a", "--scaffold_white_list", action="store", dest="scaffold_white_list", default=[],
                    type=lambda s: IdList(filename=s) if os.path.exists(s) else s.split(","),
                    help="Comma-separated list of the only scaffolds to draw. Default: all")

parser.add_argument("-b", "--scaffold_black_list", action="store", dest="scaffold_black_list", default=[],
                    type=lambda s: IdList(filename=s) if os.path.exists(s) else s.split(","),
                    help="Comma-separated list of scaffolds to skip at drawing. Default: not set")

parser.add_argument("-y", "--sort_scaffolds", action="store_true", dest="sort_scaffolds", default=False,
                    help="Order  scaffolds according to their names. Default: False")

parser.add_argument("-z", "--scaffold_ordered_list", action="store", dest="scaffold_ordered_list", default=[],
                    type=lambda s: IdList(filename=s) if os.path.exists(s) else s.split(","),
                    help="Comma-separated list of scaffolds to draw first and exactly in same order. "
                         "Scaffolds absent in this list are drawn last and in order according to vcf file . "
                         "Default: not set")
parser.add_argument("-n", "--scaffold_length_file", action="store", dest="scaffold_length_file", default=[],
                    help="File with lengths of scaffolds")
parser.add_argument("--scaffold_syn_file", action="store", dest="scaffold_syn_file",
                    help="File with scaffold id synonyms")
parser.add_argument("--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in synonym file. Default: 0")
parser.add_argument("--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in synonym file synonym. Default: 1")

"""
parser.add_argument("-q", "--figure_width", action="store", dest="figure_width", default=12, type=int,
                    help="Width of figure in inches. Default: 12")
parser.add_argument("-u", "--figure_height_scale_factor", action="store", dest="figure_height_scale_factor",
                    default=0.5, type=float,
                    help="Figure height scale factor. Figure height is calculated in inches as "
                         "int(figure_scale_factor * scaffold_number * sample_number). Default: 0.5")
"""
parser.add_argument("--masking_gff_list", action="store", dest="masking_gff_list", default=None,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of GFF files with masked regions")
parser.add_argument("--masking_threshold", action="store", dest="masking_threshold", default=0.5,
                    type=float,
                    help="Maximum gaped or masked fraction of the window. Default: 0.5")
parser.add_argument("--colormap", action="store", dest="colormap",
                    help="Matplotlib colormap to use for SNP densities. Default: not set, "
                         "colors from HapMap article are used")
parser.add_argument("--density_thresholds", action="store", dest="density_thresholds",
                    default=(0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5),
                    type=lambda s: map(float, s.split(",")),
                    help="Comma-separated list of thresholds(SNPs/kb) for SNP densities to use for window coloring. "
                         "Default: values from Hapmap article"
                         "(0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5)")

args = parser.parse_args()

variants = CollectionVCF(args.input, parsing_mode="only_coordinates")

chr_len_df = pd.DataFrame.from_csv(args.scaffold_length_file, sep='\t') if args.scaffold_length_file else variants.scaffold_length

chr_syn_dict = SynDict(filename=args.scaffold_syn_file,
                       key_index=args.syn_file_key_column,
                       value_index=args.syn_file_value_column)

if args.scaffold_syn_file:
    chr_len_df.rename(index=chr_syn_dict, inplace=True)
"""
if chr_syn_dict:
    for scaffold in raw_len_dict:
        if scaffold in chr_syn_dict:
            chr_len_df[chr_syn_dict[scaffold]] = raw_len_dict[scaffold]
    
    chr_len_df = pd.DataFrame.from_dict(chr_len_df,  orient="index")
"""
count_df = StatsVCF.count_variants_in_windows(variants, args.window_size, args.window_step,
                                              reference_scaffold_lengths=None,
                                              ignore_scaffolds_shorter_than_window=True, output_prefix=None,
                                              skip_empty_windows=False, expression=None, per_sample_output=False,
                                              scaffold_white_list=args.scaffold_white_list,
                                              scaffold_syn_dict=chr_syn_dict)

Visualization.draw_variant_window_densities(count_df, args.window_size, args.window_step, chr_len_df,
                                            args.output_prefix,
                                            figsize=(15, 10),
                                            dpi=300,
                                            colormap=args.colormap, title=args.title,
                                            extensions=args.output_formats,
                                            scaffold_order_list=args.scaffold_ordered_list)
