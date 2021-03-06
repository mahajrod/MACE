#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import pandas as pd
import argparse
from copy import deepcopy
from _collections import OrderedDict

import pandas as pd
from BCBio import GFF
from RouToolPa.Collections.General import SynDict, IdList
from RouToolPa.Parsers.VCF import CollectionVCF
from MACE.Routines import Visualization, StatsVCF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True, type=lambda s: s.split(","),
                    help="Comma_separated_list of input file with precalculated coverage in windows.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=("png", ),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg,png")

"""
parser.add_argument("-g", "--draw_gaps", action="store_true", dest="draw_gaps",
                    help="Draw gaps, ignored if reference genome is not set. Default: False")
"""
parser.add_argument("-m", "--mean_coverage_file", action="store", dest="mean_coverage_file", required=True,
                    help="File with mean coverage for all samples")
parser.add_argument("-l", "--label_list", action="store", dest="label_list", required=True, type=lambda s: s.split(","),
                    help="Comma-separated list of labels to use for samples")

parser.add_argument("--scaffold_column_name", action="store", dest="scaffold_column_name", default="scaffold",
                    help="Name of column in coverage file with scaffold ids per window. Default: scaffold")
parser.add_argument("--window_column_name", action="store", dest="window_column_name", default="window",
                    help="Name of column in coverage file with window id. Default: window")
parser.add_argument("--coverage_column_name_list", action="store", dest="coverage_column_name_list",
                    default=["median", "mean"],
                    type=lambda s: s.split(","),
                    help="Coverage file with mean/median coverage per window. Default: median,mean")
parser.add_argument("--label_column_name", action="store", dest="label_column_name", default="label",
                    help="Name of column in mean coverage file with labels of samples. Default: label")
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
                    type=lambda s: list(map(float, s.split(","))),
                    help="Comma-separated list of coverage thresholds(relative to mean/median) to use for "
                         "window coloring."
                         "Default: (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5)")
parser.add_argument("--test_colormaps", action="store_true", dest="test_colormaps",
                    help="Test colormaps. If set --colormap option will be ignored")
parser.add_argument("--absolute_coverage_values", action="store_true", dest="absolute_coverage_values",
                    help="Use absolute coverage values. Default: False")
parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float,
                    help="Adjust left border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_top", action="store", dest="subplots_adjust_top", type=float,
                    help="Adjust top border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_right", action="store", dest="subplots_adjust_right", type=float,
                    help="Adjust right border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_bottom", action="store", dest="subplots_adjust_bottom", type=float,
                    help="Adjust bottom border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--figure_width", action="store", dest="figure_width", type=float, default=15,
                    help="Width of figure in inches. Default: 15")
parser.add_argument("--figure_height_per_scaffold", action="store", dest="figure_height_per_scaffold",
                    type=float, default=0.5,
                    help="Height of figure per chromosome track. Default: 0.5")

args = parser.parse_args()

mean_coverage_df = pd.read_csv(args.mean_coverage_file, sep='\t', header=0,  index_col=0,
                               usecols=[args.label_column_name] + args.coverage_column_name_list)

chr_syn_dict = SynDict(filename=args.scaffold_syn_file,
                       key_index=args.syn_file_key_column,
                       value_index=args.syn_file_value_column)

chr_len_df = pd.read_csv(args.scaffold_length_file, sep='\t', header=None, names=("scaffold", "length"), index_col=0)

if args.scaffold_syn_file:
    chr_len_df.rename(index=chr_syn_dict, inplace=True)

coverage_df_dict = OrderedDict()

final_scaffold_set = set()
for entry, label in zip(args.input, args.label_list):
    coverage_df = pd.read_csv(entry, sep="\t", usecols=[args.scaffold_column_name,
                                                             args.window_column_name] + args.coverage_column_name_list,
                              index_col=(args.scaffold_column_name, args.window_column_name))

    scaffold_to_keep = StatsVCF.get_filtered_entry_list(coverage_df.index.get_level_values(level=0).unique().to_list(),
                                                        entry_white_list=args.scaffold_white_list)

    coverage_df = coverage_df[coverage_df.index.isin(scaffold_to_keep, level=0)]

    if args.scaffold_syn_file:
        coverage_df.rename(index=chr_syn_dict, inplace=True)

    coverage_df_dict[label] = coverage_df
    #print("AAA")
    #print(coverage_df_dict[label].index.get_level_values(0).unique().to_list())
    final_scaffold_set |= set(coverage_df_dict[label].index.get_level_values(0).unique().to_list())

#print(final_scaffold_set)
#print(chr_syn_dict)
for scaf in final_scaffold_set:
    scaf_df_list = []
    #len_df = pd.DataFrame(columns=["length"])
    len_dict = OrderedDict()
    for label in args.label_list:
        scaf_df_list.append(coverage_df_dict[label].loc[[scaf]])
        #print(scaf_df_list[-1])
        #print(label)
        scaf_df_list[-1].rename(index={scaf: label}, inplace=True)
        len_dict[label] = chr_len_df.loc[scaf]
    len_df = pd.DataFrame.from_dict(len_dict, orient="index")
    #print(len_df)
    #merged_coverage_df_dict[scaf] = pd.concat(scaf_df_list)
    #print(scaf)
    Visualization.draw_coverage_windows(pd.concat(scaf_df_list),
                                        args.window_size, args.window_step, len_df,
                                        mean_coverage_df,
                                        "{0}.{1}".format(args.output_prefix, scaf),
                                        figure_width=args.figure_width,
                                        figure_height_per_scaffold=args.figure_height_per_scaffold, dpi=300,
                                        colormap=args.colormap, title="Coverage {0}".format(scaf),
                                        extensions=args.output_formats,
                                        scaffold_order_list=args.scaffold_ordered_list,
                                        test_colormaps=args.test_colormaps,
                                        thresholds=args.coverage_thresholds,
                                        absolute_coverage_values=args.absolute_coverage_values,
                                        subplots_adjust_left=args.subplots_adjust_left,
                                        subplots_adjust_bottom=args.subplots_adjust_bottom,
                                        subplots_adjust_right=args.subplots_adjust_right,
                                        subplots_adjust_top=args.subplots_adjust_top,
                                        show_track_label=True,
                                        show_trackgroup_label=True,
                                        close_figure=True
                                        )



