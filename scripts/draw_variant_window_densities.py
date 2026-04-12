#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'


import os
import argparse
from functools import partial

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from copy import deepcopy
from RouToolPa.Collections.General import SynDict, IdList
from RouToolPa.Parsers.VCF import CollectionVCF
from RouToolPa.Parsers.BED import CollectionBED

from MACE.Routines import Visualization, StatsVCF
from MACE.Routines import Parsing
from MACE.Routines.Parsing import NewlinePreservingArgParserHelpFormatter


parser = argparse.ArgumentParser(formatter_class=NewlinePreservingArgParserHelpFormatter)
# ---- Input/Output options ----

# -------- Main Input options --------
parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input vcf file with variants or precomputed (but not colored) track.")
parser.add_argument("-t", "--input_type", action="store", dest="input_type", default="vcf",
                    help="Type of input. Allowed: 'vcf'(default), 'bedgraph'\n"
                         "\t'vcf': standard vcf format.\n"
                         "\t'bedgraph': four column bed-like format. The fourth column must contain precomputed density. Files with more columns are allowed too, but redundant columns are ignored.")

parser.add_argument("-n", "--scaffold_length_file", action="store", dest="scaffold_length_file", default=None,
                    help="File with lengths of scaffolds. Required, if input is not a VCF file.")

parser.add_argument("--centromere_bed", action="store", dest="centromere_bed", required=False,
                    type=str, help="Bed file with coordinates of centromeres")
# -------- End of Main Input options --------

# -------- Input Filtering, Renaming, Sorting and Highlighting options --------
parser.add_argument("-a", "--scaffold_whitelist", action="store", dest="scaffold_whitelist",
                    help="Comma-separated list of the only scaffolds to draw. Default: all")
parser.add_argument("--max_scaffolds", action="store", dest="max_scaffolds",
                    default=50,
                    type=int,
                    help="Maximal number of longest scaffolds from input file to show. This option works only if --scaffold_whitelist is not set. Default: 50")

parser.add_argument("-z", "--scaffold_orderlist", action="store", dest="scaffold_orderlist",
                    help="Comma-separated list of scaffolds to draw first and exactly in same order. "
                         "Scaffolds absent in this list are drawn last and in order according to vcf file . "
                         "Default: not set")

#parser.add_argument("-y", "--sort_scaffolds", action="store_true", dest="sort_scaffolds", default=False,
#                    help="Order  scaffolds according to their names. Default: False")

parser.add_argument("--scaffold_syn_file", action="store", dest="scaffold_syn_file",
                    help="File with scaffold id synonyms")
parser.add_argument("--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in synonym file. Default: 0")
parser.add_argument("--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in synonym file synonym. Default: 1")
parser.add_argument("--highlight_file", action="store", dest="highlight_file",
                    help="Tab-separated file with two columns ('scaffold' and 'color'). "
                         "Scaffold ids are ids after renaming. Must contain header.")
# -------- End of Input Filtering, Renaming and Sorting options --------

# -------- Coverage Input options --------
parser.add_argument("-c", "--coverage", action="store", dest="coverage",
                    help="File with precalculated mean/median coverage in windows")
parser.add_argument("--scaffold_column_name", action="store", dest="scaffold_column_name", default="scaffold",
                    help="Name of column in coverage file with scaffold ids per window. Default: scaffold")
parser.add_argument("--window_column_name", action="store", dest="window_column_name", default="window",
                    help="Name of column in coverage file with window id. Default: window")
parser.add_argument("--coverage_column_name", action="store", dest="coverage_column_name", default="median",
                    help="Name of column in coverage file with mean/median coverage per window. Default: median")
parser.add_argument("-m", "--mean_coverage", action="store", dest="mean_coverage",
                    type=float,
                    help="Comma-separated list of mean coverage.")
parser.add_argument("-x", "--max_coverage_threshold", action="store", dest="max_coverage_threshold", type=float, default=2.5,
                    help="Maximum coverage threshold to treat position as unmasked. Default: 2.5 * mean")
parser.add_argument("-q", "--min_coverage_threshold", action="store", dest="min_coverage_threshold", type=float, default=0.33,
                    help="Minimum coverage threshold to treat position as unmasked. Default: 0.33 * mean")
# -------- End of Coverage Input options --------

# -------- Masking Input options --------
parser.add_argument("--masking_track", action="store", dest="masking_track", default=None,
                    help="Track file with masked regions. It is 4-column bed file, the 4th column can contain "
                         "any color recognized by Matplotlib or 'masked' keyword."
                         "In case of 3-column bed file, --masking_color is used to set the masking color")
parser.add_argument("--masking_color", action="store", dest="masking_color", default='grey',
                    help="Color to use for masked windows. Default: 'grey'")

parser.add_argument("--masking_threshold", action="store", dest="masking_threshold", default=0.1,
                    type=float,
                    help="Maximum fraction of masked bases within the window. Default: 0.1")
# -------- End of Masking Input options --------

# -------- End of Input Filtering, Renaming and Sorting options --------

# -------- Output options -------
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=("png", "svg"),
                    help="Comma-separated list of formats (supported by matplotlib) for output figure. Default: svg,png")

# -------- End of Output options -------
# ---- End of Input/Output options ----

# ---- Density Track options ----
parser.add_argument("-w", "--window_size", action="store", dest="window_size", default=100000, type=int,
                    help="Size of the windows Default: 100000")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", default=None, type=int,
                    help="Step of the sliding windows. Default: window size, i.e windows are staking")
parser.add_argument("--density_multiplier", action="store", dest="density_multiplier", default=1000, type=int,
                    help="Multiplier of density. Ignored if '--prescaled_density' is set. Default: 1000, i.e. densities will be calculated per kbp")

parser.add_argument("--density_thresholds", action="store", dest="density_thresholds",
                    default=(0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5),
                    type=lambda s: list(map(float, s.split(","))),
                    help="Comma-separated list of thresholds(SNPs/kb) for SNP densities to use for window coloring. "
                         "Default: values from Hapmap article"
                         "(0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5)")
parser.add_argument("--density_thresholds_expression_type", action="store", dest="density_thresholds_expression_type",
                    default="left_open",
                    help="Type of the intervals used for thresholds. Allowed: 'left_open' (default), 'right_open' ")
parser.add_argument("--skip_top_interval", action="store_true", dest="skip_top_interval", default=False,
                    help="Skip (don't include in legend) top interval (higher than last threshold). "
                         "In such case last threshold and above values will be included in the last interval. Default: False ")
parser.add_argument("--skip_bottom_interval", action="store_true", dest="skip_bottom_interval", default=False,
                    help="Skip (don't include in legend) top interval (higher than last threshold). "
                         "In such case first threshold and below values will be included in the first interval. Default: False ")
# ---- End of Density Track  options ----

# ---- Drawing options ----

# -------- Color options --------
parser.add_argument("--only_count", action="store_true", dest="only_count", default=False,
                    help="Only count variants, do not draw them. Default: False")

parser.add_argument("--colormap", action="store", dest="colormap", default='jet',
                    help="Matplotlib colormap to use for SNP densities. Default: jet")
parser.add_argument("--custom_color_list", action="store", dest="custom_color_list", default=None,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of colors to use instead of predefined colormap. "
                         "Color names must be recognizable by matplotlib or a hex numbers preceded by #. "
                         "Number of colors in list should be equal to number of the thresholds. "
                         "If set --colormap option will be ignored. Default: not set")
parser.add_argument("--test_colormaps", action="store_true", dest="test_colormaps",
                    help="Test colormaps. If set --colormap option will be ignored")
# -------- End of Color options --------

# -------- Chromosome track options --------
parser.add_argument("--stranded", action="store_true", dest="stranded", default=False,
                    help="Stranded features and tracks. Default: False")
parser.add_argument("--rounded", action="store_true", dest="rounded", default=False,
                    help="Rounded tracks. Default: False")
parser.add_argument("--middle_break", action="store_true", dest="middle_break", default=False,
                    help="Add middle break to the track. Default: False")
parser.add_argument("--stranded_end", action="store_true", dest="stranded_end", default=False,
                    help="Stranded ends for tracks. Works only if --stranded is set. Default: False")
parser.add_argument("--feature_name", action="store", dest="feature_name", default="SNPs",
                    help="Feature name to use in legend. Default: 'SNPs'")
parser.add_argument("--fill_empty_tracks", action="store_true", dest="fill_empty_tracks", default=False,
                    help="Fill empty tracks (without features) with color set by --empty_color. Default: False")
parser.add_argument("--empty_color", action="store", dest="empty_color", default="lightgrey",
                    help="Color used to fill empty tracks. Ignored if --fill_empty_tracks is not set. "
                         "Default: 'lightgrey'")
# -------- End of Chromosome track options --------

# -------- Title options --------
parser.add_argument("-l", "--title", action="store", dest="title", default="Coverage",
                    help="Suptitle of figure. Default: 'Coverage'")
parser.add_argument("--title_fontsize", action="store", dest="title_fontsize", default=20, type=int,
                    help="Fontsize of the figure. Default: 20")
# -------- End of Title options --------

# ---- End of Drawing options ----

# -------- Label and Tick options --------
parser.add_argument("--hide_track_label", action="store_true", dest="hide_track_label", default=False,
                    help="Hide track label. Default: False")
parser.add_argument("--x_tick_type", action="store", dest="x_tick_type", default="nucleotide",
                    help="Type of xticks. Allowed: 'nucleotide' (default), 'int_number', 'float_number'")
parser.add_argument("--x_tick_fontsize", action="store", dest="x_tick_fontsize", type=int, default=None,
                    help="Fontsize of xticks. Default: matplotlib default")
# -------- End of Label and Tick options --------

# ---- Common options ----

# -------- Subplot adjustment options --------

parser.add_argument("--manual_figure_adjustment", action="store_true", dest="manual_figure_adjustment", default=False,
                    help="Adjust borders of figure manually using options below. Default: False, i.e. scaling is done automatically.")
parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float,
                    help="Adjust left border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_top", action="store", dest="subplots_adjust_top", type=float,
                    help="Adjust top border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_right", action="store", dest="subplots_adjust_right", type=float,
                    help="Adjust right border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_bottom", action="store", dest="subplots_adjust_bottom", type=float,
                    help="Adjust bottom border of subplots on the figure. Default: matplotlib defaults")
# -------- End of Subplot adjustment options --------

# -------- Figure size options --------
parser.add_argument("--figure_header_height", action="store", dest="figure_header_height",
                    type=float, default=0.0,
                    help="Height of figure header. Default: 0.0")
parser.add_argument("-u", "--figure_height_per_scaffold", action="store", dest="figure_height_per_scaffold",
                    type=float, default=0.5,
                    help="Height of figure per track. Figure height is calculated in inches as "
                         "max(1, int(figure_height_per_scaffold * scaffold_number * sample_number + figure_header_height)). Default: 0.5")
parser.add_argument("--figure_width", action="store", dest="figure_width", type=float, default=12,
                    help="Width of figure in inches. Default: 12")
# -------- End of Figure size options --------

# -------- Miscellaneous options --------
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Print additional info to stdout")
# -------- End of Miscellaneous options --------

# ---- End of Common options ----

args = parser.parse_args()
auxiliary_dict = Parsing.read_mace_auxiliary_input(len_file=args.scaffold_length_file,
                                                   whitelist_file=args.scaffold_whitelist,
                                                   max_scaffolds=args.max_scaffolds,
                                                   orderlist_file=args.scaffold_orderlist,
                                                   syn_file=args.scaffold_syn_file,
                                                   syn_file_key_column=args.syn_file_key_column,
                                                   syn_file_value_column=args.syn_file_value_column,
                                                   centromere_bed=args.centromere_bed,
                                                   highlight_bed=args.highlight_file,
                                                   legend_file=None,
                                                   vert_track_group_file=None,
                                                   hor_track_group_file=None,
                                                   hor_track_subgroup_file=None)

if args.custom_color_list is not None:
    if len(args.custom_color_list) != len(args.density_thresholds):
        raise ValueError(f"ERROR!!! Custom color list is set, but the number of colors ({len(args.custom_color_list)}) "
                         f"in the list is not equal to the number of the thresholds {len(args.density_thresholds)}!")

if args.input_type == "vcf":
    variants = CollectionVCF(args.input, parsing_mode="only_coordinates")
    if auxiliary_dict["len_df"] is None:
        auxiliary_dict["len_df"] = deepcopy(variants.scaffold_length)
        auxiliary_dict["len_df"].index.name = "scaffold"
    count_df = StatsVCF.count_variants_in_windows(variants, args.window_size, args.window_step,
                                                  reference_scaffold_lengths=auxiliary_dict["len_df"],
                                                  ignore_scaffolds_shorter_than_window=True,
                                                  output_prefix=args.output_prefix,
                                                  skip_empty_windows=False, expression=None, per_sample_output=False)

    feature_df, track_df = StatsVCF.convert_variant_count_to_feature_df(count_df,
                                                                        args.window_size,
                                                                        args.window_step)

    feature_df.to_csv("{}.features.counts".format(args.output_prefix), sep="\t", header=True, index=True)  # feature df contains a window_id (numerical) in the forth column (between end and value)
    feature_df[feature_df.columns[-1]] = feature_df[feature_df.columns[-1]] * float(args.density_multiplier) / float(args.window_size)

    feature_df.to_csv("{}.features.bed".format(args.output_prefix), sep="\t", header=True, index=True)
    track_df[track_df.columns[-1]] = track_df[track_df.columns[-1]] * float(args.density_multiplier) / float(args.window_size)

elif args.input_type in ["bedgraph"]:  # Bed format without track lines. All columns except first three are ignored. Comment lines are allowed and must start from '#'
    track_df = CollectionBED(in_file=args.input, parsing_mode="all", format="bedgraph").records
    # if args.scaffold_syn_file:
    #     track_df.rename(index=chr_syn_dict, inplace=True)

track_df = Parsing.resolve_mace_single_genome_input(track_df, auxiliary_dict)
#print(track_df)
if track_df.index.nlevels > 1:
    #drop second level of index if it was added by groupby
    track_df = track_df.groupby("scaffold").apply(lambda df: df[df["end"] <= auxiliary_dict["len_df"].loc[df.index[0], "length"]]).reset_index(level=1, drop=True)

if args.only_count:
    exit(0)

# Drawing plot
if args.custom_color_list is not None:
    cmap_list = ['custom_list']
else:
    cmap_list = Visualization.colormap_list if args.test_colormaps else [args.colormap]

for colormap in cmap_list:
    if colormap == "custom_list":
        colors = args.custom_color_list
    else:
        cmap = plt.get_cmap(colormap, len(args.density_thresholds))
        colors = [Visualization.rgb_tuple_to_hex(cmap(i)[:3]) for i in range(0, len(args.density_thresholds))]
    #print(colors)
    color_expression = partial(Visualization.color_threshold_expression,
                               thresholds=args.density_thresholds,
                               colors=colors,
                               background="white",
                               interval_type=args.density_thresholds_expression_type,
                               skip_top_interval=args.skip_top_interval,
                               skip_bottom_interval=args.skip_bottom_interval)

    track_with_colors_df = Visualization.add_color_to_track_df(track_df,
                                                               color_expression,
                                                               value_column_index=-1  # TODO fix it, add support for multiple tracks in the file
                                                               )

    track_with_colors_df.to_csv("{}.{}.track.bed".format(args.output_prefix,
                                                         colormap), sep="\t", header=True, index=True)

    Visualization.draw_features({"TR": track_with_colors_df},
                                auxiliary_dict["len_df"],
                                auxiliary_dict["orderlist_series"],
                                args.output_prefix,
                                legend=Visualization.density_legend(colors, args.density_thresholds,
                                                                    feature_name=args.feature_name,
                                                                    interval_type=args.density_thresholds_expression_type,
                                                                    skip_top_interval=args.skip_top_interval,
                                                                    skip_bottom_interval=args.skip_bottom_interval,
                                                                    masking_color=args.masking_color),
                                # query_species_color_df_dict,
                                centromere_df=auxiliary_dict["centromere_df"],
                                highlight_df=auxiliary_dict["highlight_df"],
                                figure_width=args.figure_width,
                                figure_height_per_scaffold=args.figure_height_per_scaffold,
                                figure_header_height=args.figure_header_height,
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
                                middle_break=args.middle_break,
                                stranded_end_tracks=args.stranded_end,
                                xtick_fontsize=args.x_tick_fontsize,
                                subplot_title_fontsize=args.title_fontsize,
                                subplot_title_fontweight='bold',
                                x_tick_type=args.x_tick_type,
                                autoscale_figure=False if args.manual_figure_adjustment else True,
                                fill_empty_tracks=args.fill_empty_tracks,
                                empty_color=args.empty_color,
                                )
