#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from functools import partial

import pandas as pd
import matplotlib.pyplot as plt
from RouToolPa.Collections.General import SynDict, IdList
from MACE.Routines import Visualization, StatsVCF


def rgb_tuple_to_hex(rgb_tuple):
    color_code = "#"
    for i in [0, 1, 2]:
        color_code += "{:02X}".format(int(255 * rgb_tuple[i]))

    return color_code


def read_series(s):
    return pd.read_csv(s, header=None).squeeze("columns") if os.path.exists(s) else pd.Series(s.split(","))


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with precalculated coverage in windows.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=("png", "svg"),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg, png")

parser.add_argument("-l", "--title", action="store", dest="title", default="Coverage",
                    help="Suptitle of figure. Default: 'Coverage'")
parser.add_argument("--title_fontsize", action="store", dest="title_fontsize", default=20, type=int,
                    help="Fontsize of the figure. Default: 20")
"""
parser.add_argument("-g", "--draw_gaps", action="store_true", dest="draw_gaps",
                    help="Draw gaps, ignored if reference genome is not set. Default: False")
"""
parser.add_argument("-m", "--mean_coverage_list", action="store", dest="mean_coverage_list", required=True,
                    type=lambda s: list(map(float, s.split(","))),
                    help="Comma-separated list of mean/median coverage to use")

parser.add_argument("--scaffold_column_name", action="store", dest="scaffold_column_name", default="scaffold",
                    help="Name of column in coverage file with scaffold ids per window. Default: scaffold")
parser.add_argument("--window_column_name", action="store", dest="window_column_name", default="window",
                    help="Name of column in coverage file with window id. Default: window")
parser.add_argument("--coverage_column_name_list", action="store", dest="coverage_column_name_list",
                    default=["median", "mean"],
                    type=lambda s: s.split(","),
                    help="Coverage file with mean/median coverage per window. Default: median,mean")
parser.add_argument("-w", "--window_size", action="store", dest="window_size", default=100000, type=int,
                    help="Size of the windows Default: 100000")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", default=None, type=int,
                    help="Step of the sliding windows. Default: window size, i.e windows are staking")

parser.add_argument("-a", "--scaffold_white_list", action="store", dest="scaffold_white_list", default=[],
                    type=read_series,
                    help="Comma-separated list of the only scaffolds to draw. Default: all")

parser.add_argument("-b", "--scaffold_black_list", action="store", dest="scaffold_black_list", default=[],
                    type=read_series,
                    help="Comma-separated list of scaffolds to skip at drawing. Default: not set")

parser.add_argument("-y", "--sort_scaffolds", action="store_true", dest="sort_scaffolds", default=False,
                    help="Order  scaffolds according to their names. Default: False")

parser.add_argument("-z", "--scaffold_ordered_list", action="store", dest="scaffold_ordered_list", default=[],
                    type=read_series,
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

parser.add_argument("--colormap", action="store", dest="colormap", default="jet",
                    help="Matplotlib colormap to use for coverage bins. Default: jet")
parser.add_argument("--coverage_thresholds", action="store", dest="coverage_thresholds",
                    default=(0.0, 0.25, 0.75, 1.25, 1.75, 2.25),
                    type=lambda s: list(map(float, s.split(","))),
                    help="Comma-separated list of coverage thresholds(relative to mean/median) to use for "
                         "window coloring."
                         "Default: (0.0, 0.25, 0.75, 1.25, 1.75, 2.25)")
parser.add_argument("--split_coverage_thresholds", action="store_true", dest="split_coverage_thresholds",
                    help="Use  (0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75,"
                         "0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625,"
                         "1.75, 1.875, 2.0, 2.125, 2.25) as thresholds instead  of default ones."
                         "Doesn't work if --coverage_thresholds is set")
parser.add_argument("--test_colormaps", action="store_true", dest="test_colormaps",
                    help="Test colormaps. If set --colormap option will be ignored")

parser.add_argument("--hide_track_label", action="store_true", dest="hide_track_label", default=False,
                    help="Hide track label. Default: False")

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
parser.add_argument("--x_tick_fontsize", action="store", dest="x_tick_fontsize", type=int, default=None,
                    help="Fontsize of xticks. Default: matplotlib default")
parser.add_argument("--stranded", action="store_true", dest="stranded", default=False,
                    help="Stranded features and tracks. Default: False")
parser.add_argument("--rounded", action="store_true", dest="rounded", default=False,
                    help="Rounded tracks. Default: False")
parser.add_argument("--stranded_end", action="store_true", dest="stranded_end", default=False,
                    help="Stranded ends for tracks. Works only if --stranded is set. Default: False")
parser.add_argument("--centromere_bed", action="store", dest="centromere_bed", required=False,
                    type=str,
                    help="Bed file with coordinates of centromeres")
parser.add_argument("--highlight_file", action="store", dest="highlight_file",
                    type=lambda s: pd.read_csv(s, header=0, index_col=0, sep="\t"),
                    help="Tab-separated file with two columns ('scaffold' and 'color'). "
                         "Scaffold ids are ids after renaming"
                         "Must contain header.")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Print additional info to stdout")

args = parser.parse_args()

if args.window_step is None:
    args.window_step = args.window_size

args.scaffold_ordered_list = args.scaffold_ordered_list[::-1]

chr_syn_dict = SynDict(filename=args.scaffold_syn_file,
                       key_index=args.syn_file_key_column,
                       value_index=args.syn_file_value_column)

if args.centromere_bed:
    centromere_df = pd.read_csv(args.centromere_bed,
                                usecols=(0, 1, 2),
                                index_col=0,
                                header=None,
                                sep="\t", names=["scaffold_id", "start", "end"])
    centromere_df.rename(index=chr_syn_dict, inplace=True)
else:
    centromere_df = None


coverage_df = pd.read_csv(args.input, sep="\t", usecols=[args.scaffold_column_name,
                                                         args.window_column_name] + args.coverage_column_name_list,
                          index_col=(args.scaffold_column_name, args.window_column_name),)

coverage_df.index.set_levels(list(map(str, coverage_df.index.levels[0])), level=0, inplace=True)
if args.verbose:
    print("Coverage df (raw)")
    print(coverage_df)
    print("Coverage df index")
    print(coverage_df.index.get_level_values(level=0))
    print("Whitelist")
    print(args.scaffold_white_list)
scaffold_to_keep = StatsVCF.get_filtered_entry_list(coverage_df.index.get_level_values(level=0).unique().to_list(),
                                                    entry_white_list=args.scaffold_white_list)

coverage_df = coverage_df[coverage_df.index.isin(scaffold_to_keep, level=0)]
chr_len_df = pd.read_csv(args.scaffold_length_file, sep='\t', header=None, names=("scaffold", "length"), index_col=0,
                         converters={args.scaffold_column_name: str,})
chr_len_df.index = pd.Index(list(map(str, chr_len_df.index)))

if args.scaffold_syn_file:
    coverage_df.rename(index=chr_syn_dict, inplace=True)
    chr_len_df.rename(index=chr_syn_dict, inplace=True)

average_coverage_dict = dict(zip(args.coverage_column_name_list, args.mean_coverage_list))

if args.verbose:
    print(chr_syn_dict)
    print(scaffold_to_keep)
    print(coverage_df)
    print(chr_syn_dict)
    print(coverage_df)

#print(coverage_df)

args.coverage_thresholds = (0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75,
                            0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625,
                            1.75, 1.875, 2.0, 2.125, 2.25) if args.split_coverage_thresholds else args.coverage_thresholds

track_df_dict = {}
cmap = plt.get_cmap(args.colormap, len(args.coverage_thresholds))
colors = [rgb_tuple_to_hex(cmap(i)[:3]) for i in range(0, len(args.coverage_thresholds))]


for mean_coverage, track_label in zip(args.mean_coverage_list, args.coverage_column_name_list):
    feature_df, track_df = StatsVCF.convert_variant_count_to_feature_df(coverage_df,
                                                                        args.window_size,
                                                                        args.window_step,
                                                                        value_column=track_label)
    track_df[track_label] /= mean_coverage
    color_expression = partial(Visualization.color_threshold_expression,
                               thresholds=args.coverage_thresholds,
                               colors=colors,
                               background="white")

    track_with_colors_df = Visualization.add_color_to_track_df(track_df,
                                                               color_expression,
                                                               value_column_index=-1 # TODO fix it, add support for multiple tracks in the file
                                                               )
    #print(track_with_colors_df)
    track_df_dict[track_label] = track_with_colors_df

Visualization.draw_features(track_df_dict,
                            chr_len_df,
                            args.scaffold_ordered_list,
                            args.output_prefix,
                            legend=Visualization.coverage_legend(colormap=args.colormap, thresholds=args.coverage_thresholds), #Visualization.density_legend(colors, args.coverage_thresholds),
                            # query_species_color_df_dict,
                            centromere_df=centromere_df,
                            highlight_df=args.highlight_file,
                            figure_width=args.figure_width,
                            figure_height_per_scaffold=args.figure_height_per_scaffold,
                            dpi=300,
                            default_color="red",
                            title=args.title,
                            extensions=args.output_formats,
                            feature_start_column_id="start",
                            feature_end_column_id="end",
                            feature_color_column_id="color",
                            feature_length_column_id="length",
                            subplots_adjust_left=args.subplots_adjust_left,
                            subplots_adjust_bottom=args.subplots_adjust_bottom,
                            subplots_adjust_right=args.subplots_adjust_right,
                            subplots_adjust_top=args.subplots_adjust_top,
                            show_track_label=not args.hide_track_label,
                            show_trackgroup_label=True,
                            close_figure=True,
                            subplot_scale=False,
                            track_group_scale=False,
                            track_group_distance=2,
                            xmax_multiplier=1.3, ymax_multiplier=1.00,
                            stranded_tracks=args.stranded,
                            rounded_tracks=args.rounded,
                            stranded_end_tracks=args.stranded_end,
                            xtick_fontsize=args.x_tick_fontsize,
                            subplot_title_fontsize=args.title_fontsize,
                            subplot_title_fontweight='bold'
                            )
"""    
Visualization.draw_coverage_windows(coverage_df, args.window_size, args.window_step, chr_len_df,
                                    average_coverage_dict,
                                    args.output_prefix,
                                    figure_width=args.figure_width,
                                    figure_height_per_scaffold=args.figure_height_per_scaffold, dpi=300,
                                    colormap=args.colormap, title=args.title,
                                    extensions=args.output_formats,
                                    scaffold_order_list=args.scaffold_ordered_list,
                                    test_colormaps=args.test_colormaps,
                                    thresholds=(0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75,
                                                0.875, 1.0, 1.125, 1.25, 1.375, 1.5, 1.625,
                                                1.75, 1.875, 2.0, 2.125, 2.25) if args.split_coverage_thresholds else args.coverage_thresholds,
                                    absolute_coverage_values=args.absolute_coverage_values,
                                    subplots_adjust_left=args.subplots_adjust_left,
                                    subplots_adjust_bottom=args.subplots_adjust_bottom,
                                    subplots_adjust_right=args.subplots_adjust_right,
                                    subplots_adjust_top=args.subplots_adjust_top,
                                    show_track_label=not args.hide_track_label,
                                    show_trackgroup_label=True
                                    )

"""

