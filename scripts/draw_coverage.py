#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import pandas as pd
import argparse
from copy import deepcopy

import pandas as pd
from BCBio import GFF
from RouToolPa.Collections.General import SynDict, IdList
from RouToolPa.Parsers.VCF import CollectionVCF
from MACE.Routines import Visualization, StatsVCF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with precalculated coverage in windows.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=("png", ),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg,png")

parser.add_argument("-l", "--title", action="store", dest="title", default="Coverage",
                    help="Suptitle of figure. Default: 'Coverage'")
"""
parser.add_argument("-g", "--draw_gaps", action="store_true", dest="draw_gaps",
                    help="Draw gaps, ignored if reference genome is not set. Default: False")
"""
parser.add_argument("-m", "--mean_coverage", action="store", dest="mean_coverage", required=True, type=float,
                    help="Mean/median coverage to use")

parser.add_argument("--scaffold_column_name", action="store", dest="scaffold_column_name", default="scaffold",
                    help="Name of column in coverage file with scaffold ids per window. Default: scaffold")
parser.add_argument("--window_column_name", action="store", dest="window_column_name", default="window",
                    help="Name of column in coverage file with window id. Default: window")
parser.add_argument("--coverage_column_name", action="store", dest="coverage_column_name", default="median",
                    help="Name of column in coverage file with mean/median coverage per window. Default: median")

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
parser.add_argument("-n", "--scaffold_length_file", action="store", dest="scaffold_length_file", required=True,
                    help="File with lengths of scaffolds")
parser.add_argument("--scaffold_syn_file", action="store", dest="scaffold_syn_file",
                    help="File with scaffold id synonyms")
parser.add_argument("--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in synonym file. Default: 0")
parser.add_argument("--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in synonym file synonym. Default: 1")

parser.add_argument("--colormap", action="store", dest="colormap",
                    help="Matplotlib colormap to use for SNP densities. Default: not set, "
                         "colors from HapMap article are used")
parser.add_argument("--coverage_thresholds", action="store", dest="coverage_thresholds",
                    default=(0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5),
                    type=lambda s: map(float, s.split(",")),
                    help="Comma-separated list of coverage thresholds(relative to mean/median) to use for "
                         "window coloring."
                         "(0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5)")
parser.add_argument("--test_colormaps", action="store_true", dest="test_colormaps",
                    help="Test colormaps. If set --colormap option will be ignored")

args = parser.parse_args()

chr_syn_dict = SynDict(filename=args.scaffold_syn_file,
                       key_index=args.syn_file_key_column,
                       value_index=args.syn_file_value_column)


coverage_df = pd.read_csv(args.input, sep="\t", usecols=(args.scaffold_column_name,
                                                         args.window_column_name,
                                                         args.coverage_column_name),
                          index_col=(args.scaffold_column_name, args.window_column_name))

scaffold_to_keep = StatsVCF.get_filtered_entry_list(coverage_df.index.get_level_values(level=0).unique().to_list(),
                                                    entry_white_list=args.scaffold_white_list)
#print(scaffold_to_keep)
coverage_df = coverage_df[coverage_df.index.isin(scaffold_to_keep, level=0)]

chr_len_df = pd.read_csv(args.scaffold_length_file, sep='\t', header=None, names=("scaffold", "length"), index_col=0)


#print(coverage_df)

if args.scaffold_syn_file:
    coverage_df.rename(index=chr_syn_dict, inplace=True)
    chr_len_df.rename(index=chr_syn_dict, inplace=True)

#print(chr_syn_dict)
#print(coverage_df)
#print(coverage_df.index)
#print(chr_len_df)
Visualization.draw_coverage_windows(coverage_df, args.window_size, args.window_step, chr_len_df,
                                    args.mean_coverage,
                                    args.output_prefix,
                                    figure_width=15, figure_height_per_scaffold=0.5, dpi=300,
                                    colormap=args.colormap, title=args.title,
                                    extensions=args.output_formats,
                                    scaffold_order_list=args.scaffold_ordered_list,
                                    test_colormaps=args.test_colormaps,
                                    thresholds=args.coverage_thresholds)



