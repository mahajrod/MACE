#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import argparse
from functools import partial

import pandas as pd

import matplotlib.pyplot as plt
from copy import deepcopy
from RouToolPa.Collections.General import SynDict, IdList
from RouToolPa.Parsers.VCF import CollectionVCF
from MACE.Routines import Visualization, StatsVCF


def read_series(s):
    return pd.read_csv(s, header=None).squeeze("columns") if os.path.exists(s) else pd.Series(s.split(","))


def rgb_tuple_to_hex(rgb_tuple):
    color_code = "#"
    for i in [0, 1, 2]:
        color_code += "{:02X}".format(int(255 * rgb_tuple[i]))

    return color_code


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input vcf file with variants or precomputed track.")
parser.add_argument("-t", "--input_type", action="store", dest="input_type", default="vcf",
                    help="Type of input. Allowed: 'vcf'(default), 'bedgraph'")
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
                    default=("png", "svg"),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg,png")

parser.add_argument("-l", "--title", action="store", dest="title", default="Variant density",
                    help="Suptitle of figure. Default: 'Variant density'")
parser.add_argument("--title_fontsize", action="store", dest="title_fontsize", default=20, type=int,
                    help="Fontsize of the figure. Default: 20")

"""
parser.add_argument("-g", "--draw_gaps", action="store_true", dest="draw_gaps",
                    help="Draw gaps, ignored if reference genome is not set. Default: False")
"""
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
parser.add_argument("-q", "--min_coverage_threshold", action="store", dest="min_coverage_threshold", type=float, default=0.5,
                    help="Minimum coverage threshold to treat position as unmasked. Default: 0.5 * mean")

parser.add_argument("-w", "--window_size", action="store", dest="window_size", default=100000, type=int,
                    help="Size of the windows Default: 100000")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", default=None, type=int,
                    help="Step of the sliding windows. Default: window size, i.e windows are staking")
parser.add_argument("--density_multiplier", action="store", dest="density_multiplier", default=1000, type=int,
                    help="Multiplier of density. Default: 1000, i.e. densities will be calculated per kbp")
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


parser.add_argument("--figure_width", action="store", dest="figure_width", default=12, type=int,
                    help="Width of figure in inches. Default: 12")
parser.add_argument("-u", "--figure_height_per_scaffold", action="store", dest="figure_height_per_scaffold",
                    default=0.5, type=float,
                    help="Figure height per scaffold. Figure height is calculated in inches as "
                         "int(figure_height_per_scaffold * scaffold_number * sample_number). Default: 0.5")
parser.add_argument("--figure_header_height", action="store", dest="figure_header_height",
                    type=float, default=0,
                    help="Height of figure header. Default: 0")
parser.add_argument("--masking_gff_list", action="store", dest="masking_gff_list", default=None,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of GFF files with masked regions")
parser.add_argument("--masking_threshold", action="store", dest="masking_threshold", default=0.5,
                    type=float,
                    help="Maximum gaped or masked fraction of the window. Default: 0.5")
parser.add_argument("--colormap", action="store", dest="colormap", default='jet',
                    help="Matplotlib colormap to use for SNP densities. Default: jet")
parser.add_argument("--density_thresholds", action="store", dest="density_thresholds",
                    default=(0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5),
                    type=lambda s: list(map(float, s.split(","))),
                    help="Comma-separated list of thresholds(SNPs/kb) for SNP densities to use for window coloring. "
                         "Default: values from Hapmap article"
                         "(0.0, 0.1, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5)")
parser.add_argument("--test_colormaps", action="store_true", dest="test_colormaps",
                    help="Test colormaps. If set --colormap option will be ignored")
parser.add_argument("--hide_track_label", action="store_true", dest="hide_track_label", default=False,
                    help="Hide track label. Default: False")
parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float,
                    help="Adjust left border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_top", action="store", dest="subplots_adjust_top", type=float,
                    help="Adjust top border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_right", action="store", dest="subplots_adjust_right", type=float,
                    help="Adjust right border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_bottom", action="store", dest="subplots_adjust_bottom", type=float,
                    help="Adjust bottom border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--only_count", action="store_true", dest="only_count", default=False,
                    help="Only count variants, do not draw them. Default: False")
parser.add_argument("--x_tick_fontsize", action="store", dest="x_tick_fontsize", type=int, default=None,
                    help="Fontsize of xticks. Default: matplotlib default")
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
parser.add_argument("--centromere_bed", action="store", dest="centromere_bed", required=False,
                    type=str,
                    help="Bed file with coordinates of centromeres")
parser.add_argument("--highlight_file", action="store", dest="highlight_file",
                    type=lambda s: pd.read_csv(s, header=0, index_col=0, sep="\t"),
                    help="Tab-separated file with two columns ('scaffold' and 'color'). "
                         "Scaffold ids are ids after renaming"
                         "Must contain header.")
args = parser.parse_args()

args.scaffold_ordered_list = args.scaffold_ordered_list[::-1]

if isinstance(args.scaffold_ordered_list, list):
    if not args.scaffold_ordered_list:
        args.scaffold_ordered_list = args.scaffold_white_list
else:
    if args.scaffold_ordered_list.empty:
        args.scaffold_ordered_list = args.scaffold_white_list

variants = CollectionVCF(args.input, parsing_mode="only_coordinates")

chr_len_df = pd.read_csv(args.scaffold_length_file, sep='\t', header=None, index_col=0) if args.scaffold_length_file else deepcopy(variants.scaffold_length)
chr_len_df.index = pd.Index(list(map(str, chr_len_df.index)))
chr_len_df.index.name = "scaffold"
chr_len_df.columns = ["length"]

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
#print(chr_syn_dict)

if args.input_type == "vcf":
    count_df = StatsVCF.count_variants_in_windows(variants, args.window_size, args.window_step,
                                                  reference_scaffold_lengths=chr_len_df,
                                                  ignore_scaffolds_shorter_than_window=True,
                                                  output_prefix=args.output_prefix,
                                                  skip_empty_windows=False, expression=None, per_sample_output=False,
                                                  scaffold_white_list=args.scaffold_white_list,
                                                  scaffold_syn_dict=chr_syn_dict)
    feature_df, track_df = StatsVCF.convert_variant_count_to_feature_df(count_df,
                                                                    args.window_size,
                                                                    args.window_step)
    feature_df.to_csv("{}.features.counts".format(args.output_prefix), sep="\t", header=True, index=True)
    feature_df[feature_df.columns[-1]] = feature_df[feature_df.columns[-1]] * float(args.density_multiplier) / float(args.window_size)

    feature_df.to_csv("{}.features.bed".format(args.output_prefix), sep="\t", header=True, index=True)

elif args.input_type == "bedgraph":
    track_df = pd.read_csv(args.input, sep="\t", names=["scaffold", "start", "end", "value"],
                           header=None, index_col=0, na_values=".")
    track_df["value"] = track_df["value"].astype(float)

if args.scaffold_syn_file:
    chr_len_df.rename(index=chr_syn_dict, inplace=True)

track_df[track_df.columns[-1]] = track_df[track_df.columns[-1]] * float(args.density_multiplier) / float(args.window_size)


# TODO: rewrite application of masking
"""
if args.coverage:
    masking_df = pd.read_csv(args.coverage, usecols=(args.scaffold_column_name,
                                                     args.window_column_name,
                                                     args.coverage_column_name),
                             index_col=(args.scaffold_column_name, args.window_column_name),
                             sep="\t")
    scaffold_to_keep = StatsVCF.get_filtered_entry_list(
        masking_df.index.get_level_values(level=0).unique().to_list(),
        entry_white_list=args.scaffold_white_list)
    masking_df = masking_df[masking_df.index.isin(scaffold_to_keep, level=0)]
    # print(scaffold_to_keep)
    if chr_syn_dict:
        masking_df.rename(index=chr_syn_dict, inplace=True)
    # print("aaaaaaaa")
    # print(masking_df)

    min_threshold = args.mean_coverage * args.min_coverage_threshold
    max_threshold = args.mean_coverage * args.max_coverage_threshold
    count_df["masked"] = (masking_df[args.coverage_column_name] < min_threshold) | (
            masking_df[args.coverage_column_name] > max_threshold)
    # print(masking_df)
    # print(count_df)
    count_df.to_csv("%s.variant_counts.with_masking.tsv" % args.output_prefix, sep='\t', header=True, index=True)
"""
if not args.only_count:

    for colormap in (Visualization.colormap_list if args.test_colormaps else [args.colormap]):
        #title = (args.title + " (colormap {})".format(colormap)) if args.title else "Colormap {}".format(colormap)

        cmap = plt.get_cmap(colormap, len(args.density_thresholds))
        colors = [rgb_tuple_to_hex(cmap(i)[:3]) for i in range(0, len(args.density_thresholds))]

        color_expression = partial(Visualization.color_threshold_expression,
                                   thresholds=args.density_thresholds,
                                   colors=colors,
                                   background="white")

        track_with_colors_df = Visualization.add_color_to_track_df(track_df,
                                                                   color_expression,
                                                                   value_column_index=-1 # TODO fix it, add support for multiple tracks in the file
                                                                   )

        track_with_colors_df.to_csv("{}.{}.track.bed".format(args.output_prefix,
                                                               colormap), sep="\t", header=True, index=True)
        #print(feature_with_colors_df)
        #print(args.scaffold_ordered_list)
        Visualization.draw_features({"TR": track_with_colors_df},
                                    chr_len_df,
                                    args.scaffold_ordered_list,
                                    args.output_prefix,
                                    legend=Visualization.density_legend(colors, args.density_thresholds,
                                                                        feature_name=args.feature_name),
                                    # query_species_color_df_dict,
                                    centromere_df=centromere_df,
                                    highlight_df=args.highlight_file,
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
                                    subplot_title_fontweight='bold'
                                    )
        """
        Visualization.draw_variant_window_densities(count_df, args.window_size, args.window_step, chr_len_df,
                                                    args.output_prefix,
                                                    figure_width=15,
                                                    figure_height_per_scaffold=0.5,
                                                    dpi=300,
                                                    show_track_label=not args.hide_track_label,
                                                    colormap=args.colormap, title=args.title,
                                                    extensions=args.output_formats,
                                                    scaffold_order_list=args.scaffold_ordered_list,
                                                    test_colormaps=args.test_colormaps,
                                                    thresholds=args.density_thresholds,
                                                    masking=True if args.coverage else False,
                                                    subplots_adjust_left=args.subplots_adjust_left,
                                                    subplots_adjust_bottom=args.subplots_adjust_bottom,
                                                    subplots_adjust_right=args.subplots_adjust_right,
                                                    subplots_adjust_top=args.subplots_adjust_top,
                                                    show_trackgroup_label=True)
        """