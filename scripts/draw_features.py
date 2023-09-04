#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os

import argparse

import pandas as pd
import numpy as np
from RouToolPa.Parsers.STR import CollectionSTR
from RouToolPa.Parsers.GFF import CollectionGFF
from RouToolPa.Parsers.BLAST import CollectionBLAST
from RouToolPa.Parsers.BED import CollectionBED
from RouToolPa.Collections.General import SynDict, IdList
from MACE.Routines import Visualization, StatsVCF


def read_series(s):
    return pd.read_csv(s, header=None).squeeze("columns") if os.path.exists(s) else pd.Series(s.split(","))


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with selected features")
parser.add_argument("-t", "--input_type", action="store", dest="input_type", default="str",
                    help="Type of input file. Allowed: str (default), gff, blast, bed, bedgraph")
parser.add_argument("-g", "--legend", action="store", dest="legend",
                    help="File with legend for feature colors containing two columns with color and legend text")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=("png", "svg"),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg,png")

parser.add_argument("-l", "--title", action="store", dest="title", default="Coverage",
                    help="Suptitle of figure. Default: 'Coverage'")

parser.add_argument("--scaffold_column_name", action="store", dest="scaffold_column_name",
                    help="Name of column in feature file with scaffold ids . Default: dependent on format of the file")
parser.add_argument("--start_column_name", action="store", dest="start_column_name",
                    help="Name of column in feature file with starts. Default: dependent on format of the file")
parser.add_argument("--end_column_name", action="store", dest="end_column_name",
                    help="Name of column in feature file with ends. Default: dependent on format of the file")
parser.add_argument("--color_column_name", action="store", dest="color_column_name",
                    help="Name of column in feature file with color. Default: not set")
parser.add_argument("--default_color", action="store", dest="default_color", default="red",
                    help="Default color used for all features if color column is not set. Default: red")
parser.add_argument("-a", "--scaffold_white_list", action="store", dest="scaffold_white_list",
                    default=pd.Series(dtype=str),
                    type=read_series,
                    help="Comma-separated list of the only scaffolds to draw. Default: all")

parser.add_argument("-b", "--scaffold_black_list", action="store", dest="scaffold_black_list",
                    default=pd.Series(dtype=str),
                    type=read_series,
                    help="Comma-separated list of scaffolds to skip at drawing. Default: not set")

parser.add_argument("-y", "--sort_scaffolds", action="store_true", dest="sort_scaffolds", default=False,
                    help="Order  scaffolds according to their names. Default: False")

parser.add_argument("-z", "--scaffold_ordered_list", action="store", dest="scaffold_ordered_list",
                    default=pd.Series(dtype=str),
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
                    help="Matplotlib colormap to use for SNP densities. Default: jet")


parser.add_argument("--hide_track_label", action="store_true", dest="hide_track_label", default=False,
                    help="Hide track label. Default: False")
parser.add_argument("--x_tick_type", action="store", dest="x_tick_type", default="nucleotide",
                    help="Type of xticks. Allowed: 'nucleotide' (default), 'int_number', 'float_number'")

parser.add_argument("--feature_shape", action="store", dest="feature_shape", default="rectangle",
                    help="Shape of features. Allowed: rectangle(default), circle, ellipse")

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
parser.add_argument("--figure_header_height", action="store", dest="figure_header_height",
                    type=float, default=0.0,
                    help="Height of figure header. Default: 0.0")
parser.add_argument("--figure_height_per_scaffold", action="store", dest="figure_height_per_scaffold",
                    type=float, default=0.5,
                    help="Height of figure per chromosome track. Default: 0.5")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Print additional info to stdout")

parser.add_argument("--subplot_scale", action="store_true", dest="subplot_scale",
                    help="Scale feature x size by subplot x/y ratio. Default: off")
parser.add_argument("--track_group_scale", action="store_true", dest="track_group_scale",
                    help="Scale feature x size by track_group x/y ratio. Default: off")
parser.add_argument("--x_tick_fontsize", action="store", dest="x_tick_fontsize", type=int, default=None,
                    help="Fontsize of xticks. Default: matplotlib default")
parser.add_argument("--stranded", action="store_true", dest="stranded", default=False,
                    help="Stranded features and tracks. Default: False")
parser.add_argument("--rounded", action="store_true", dest="rounded", default=False,
                    help="Rounded tracks. Default: False")
parser.add_argument("--stranded_end", action="store_true", dest="stranded_end", default=False,
                    help="Stranded ends for tracks. Works only if --stranded is set. Default: False")
parser.add_argument("--fill_empty_tracks", action="store_true", dest="fill_empty_tracks", default=False,
                    help="Fill empty tracks (without features) with color set by --empty_color. Default: False")
parser.add_argument("--empty_color", action="store", dest="empty_color", default="lightgrey",
                    help="Color used to fill empty tracks. Ignored if --fill_empty_tracks is not set. "
                         "Default: 'lightgrey'")
parser.add_argument("--centromere_bed", action="store", dest="centromere_bed", required=False,
                    type=str,
                    help="Bed file with coordinates of centromeres")
parser.add_argument("--highlight_file", action="store", dest="highlight_file",
                    type=lambda s: pd.read_csv(s, header=0, index_col=0, sep="\t"),
                    help="Tab-separated file with two columns ('scaffold' and 'color'). "
                         "Scaffold ids are ids after renaming"
                         "Must contain header.")
parser.add_argument("--title_fontsize", action="store", dest="title_fontsize", default=20, type=int,
                    help="Fontsize of the figure. Default: 20")

args = parser.parse_args()

args.scaffold_ordered_list = args.scaffold_ordered_list[::-1]

if isinstance(args.scaffold_ordered_list, list):
    if not args.scaffold_ordered_list:
        args.scaffold_ordered_list = args.scaffold_white_list
else:
    if args.scaffold_ordered_list.empty:
        args.scaffold_ordered_list = args.scaffold_white_list

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
    print(centromere_df)
else:
    centromere_df = None
try:
    if args.input_type == "str":
        feature_df = CollectionSTR(in_file=args.input, records=None, format="filtered_str", parsing_mode="all",
                                   black_list=(), white_list=(),
                                   syn_dict=chr_syn_dict)

        feature_df.records.set_index("scaffold_id", inplace=True)

        feature_start_column_id = args.start_column_name if args.start_column_name else "start"
        feature_end_column_id = args.end_column_name if args.end_column_name else "end"

    elif args.input_type == "gff":
        feature_df = CollectionGFF(in_file=args.input, parsing_mode="only_coordinates")

        feature_start_column_id = args.start_column_name if args.start_column_name else "start"
        feature_end_column_id = args.end_column_name if args.end_column_name else "end"

    elif args.input_type == "bed":
        feature_df = CollectionBED(in_file=args.input, parsing_mode="coordinates_only", format="bed")

        feature_start_column_id = "start"
        feature_end_column_id = "end"

    elif args.input_type == "bedgraph":
        feature_df = CollectionBED(in_file=args.input, parsing_mode="all", format="bed")
        feature_df.records.columns = pd.Index(["start", "end", "value"])
        feature_df.records.index.name = "scaffold"
        feature_start_column_id = "start"
        feature_end_column_id = "end"

    elif args.input_type == "blast":
        feature_df = CollectionBLAST(in_file=args.input, parsing_mode="complete")
        feature_df.records.reset_index(level="query_id", inplace=True)
        feature_start_column_id = args.start_column_name if args.start_column_name else "target_start"
        feature_end_column_id = args.end_column_name if args.end_column_name else "target_end"

    else:
        raise ValueError("ERROR!!! Unrecognized input type ({}). ".format(args.input_type))
except pd.errors.EmptyDataError:
    print("Empty input file. Silent exit.")   # try-except added to handle case when input file is empty without raising exception. For use in snakemake
    exit(0)

legend_df = pd.read_csv(args.legend, header=None, index_col=0, sep="\t") if args.legend else None


scaffold_to_keep = StatsVCF.get_filtered_entry_list(feature_df.records.index.get_level_values(level=0).unique().to_list(),
                                                    entry_white_list=args.scaffold_white_list)

#print(scaffold_to_keep)
# remove redundant scaffolds
feature_df.records = feature_df.records[feature_df.records.index.isin(scaffold_to_keep)]
args.scaffold_ordered_list = args.scaffold_ordered_list[args.scaffold_ordered_list.isin(pd.Series(args.scaffold_white_list).replace(chr_syn_dict))]
#print(scaffold_to_keep)
#print(pd.Series(scaffold_to_keep).replace(chr_syn_dict))
chr_len_df = pd.read_csv(args.scaffold_length_file, sep='\t', header=None, names=("scaffold", "length"), index_col=0)
chr_len_df.index = pd.Index(list(map(str, chr_len_df.index)))

if args.scaffold_syn_file:
    chr_len_df.rename(index=chr_syn_dict, inplace=True)
    feature_df.records.rename(index=chr_syn_dict, inplace=True)
if args.verbose:
    print(chr_syn_dict)
    print(feature_df.records)
#print(feature_df.records.columns)
#print(feature_df.records)
#print(chr_len_df)

Visualization.draw_features({"features": feature_df}, chr_len_df,
                            args.scaffold_ordered_list,
                            args.output_prefix,
                            legend=Visualization.feature_legend(legend_df, colormap=args.colormap),
                            #legend_df=legend_df,
                            centromere_df=centromere_df,
                            highlight_df=args.highlight_file,
                            figure_width=args.figure_width,
                            figure_height_per_scaffold=0.5,
                            figure_header_height=args.figure_header_height,
                            dpi=300,
                            #colormap=None, thresholds=None, colors=None, background=None,
                            default_color=args.default_color,
                            title=args.title,
                            extensions=args.output_formats,
                            feature_shape=args.feature_shape,
                            feature_start_column_id=feature_start_column_id,
                            feature_end_column_id=feature_end_column_id,
                            feature_color_column_id=args.color_column_name,
                            feature_length_column_id="length",
                            subplots_adjust_left=args.subplots_adjust_left,
                            subplots_adjust_bottom=args.subplots_adjust_bottom,
                            subplots_adjust_right=args.subplots_adjust_right,
                            subplots_adjust_top=args.subplots_adjust_top,
                            show_track_label=not args.hide_track_label,
                            show_trackgroup_label=True,
                            subplot_scale=args.subplot_scale,
                            track_group_scale=args.track_group_scale,
                            stranded_tracks=args.stranded,
                            rounded_tracks=args.rounded,
                            stranded_end_tracks=args.stranded_end,
                            fill_empty_tracks=args.fill_empty_tracks,
                            empty_color=args.empty_color,
                            xtick_fontsize=args.x_tick_fontsize,
                            subplot_title_fontsize=args.title_fontsize,
                            subplot_title_fontweight='bold',
                            x_tick_type=args.x_tick_type,
                            )



