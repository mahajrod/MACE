#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os

import argparse
from copy import deepcopy
from collections import OrderedDict
import pandas as pd
import numpy as np

from distinctipy import distinctipy
from functools import partial

import xlsxwriter as xlsx
from xlsxwriter.utility import xl_rowcol_to_cell

from RouToolPa.Collections.General import SynDict, IdList
from MACE.Routines import Visualization, StatsVCF

from RouToolPa.Parsers.PSL import CollectionPSL
from RouToolPa.Parsers.BED import CollectionBED


def split_comma_separated_list(string):
    return string.split(",")


def rgb_tuple_to_hex(rgb_tuple):
    color_code = "#"
    for i in [0, 1, 2]:
        color_code += "{:02X}".format(int(255 * rgb_tuple[i]))

    return color_code


def dist_to_next_seg_from_same_chr(df):
    temp = deepcopy(df)
    temp["dist_upstream"] = df["start"].astype("Int64") - df["end"].astype("Int64").shift(1)
    temp["dist_downstream"] = df["start"].astype("Int64").shift(-1) - df["end"].astype("Int64")

    temp["dist_end_start"] = df["start"].astype("Int64") - df["end"].astype("Int64").shift(1)
    temp["dist_end_end"] = df["end"].astype("Int64") - df["end"].astype("Int64").shift(1)
    temp["embedded"] = temp["dist_end_end"] <= 0
    return temp


def filter_isolated_short_blocks(bed_df, min_block_len=1000000, max_dist_between_short_blocks=3000000):
    tmp_df = bed_df.reset_index(drop=False).set_index(["scaffold", "query"]).sort_values(by=["scaffold",
                                                                                             "start",
                                                                                             "end",
                                                                                             "query",
                                                                                             "query_start",
                                                                                             "query_end"], axis=0)
    tmp_df["target_len"] = tmp_df["end"] - tmp_df["start"]
    tmp_df["query_len"] = tmp_df["query_end"] - tmp_df["query_start"]

    tmp_df = tmp_df.groupby(["scaffold",
                             "query"], sort=False).apply(dist_to_next_seg_from_same_chr).sort_values(by=["scaffold",
                                                                                                         "start",
                                                                                                         "end",
                                                                                                         "query",
                                                                                                         "query_start",
                                                                                                         "query_end"],
                                                                                                     axis=0)  # .reset_index(level=(0,1), drop=True)
    tmp_df["embedded"].fillna(0, inplace=True)
    # recalculate distances after removal of embedded blocks
    tmp_df = tmp_df[tmp_df["embedded"] <= 0]
    tmp_df = tmp_df.groupby(["scaffold",
                             "query"], sort=False).apply(dist_to_next_seg_from_same_chr).sort_values(by=["scaffold",
                                                                                                         "start",
                                                                                                         "end",
                                                                                                         "query",
                                                                                                         "query_start",
                                                                                                         "query_end"],
                                                                                                     axis=0)
    tmp_df = tmp_df[~((tmp_df["target_len"] < min_block_len) & (((tmp_df["dist_upstream"] > max_dist_between_short_blocks) | tmp_df["dist_upstream"].isna()) & (
            (tmp_df["dist_downstream"] > max_dist_between_short_blocks) | (tmp_df["dist_downstream"].isna()))))]

    return tmp_df


def merge_adjacent_blocks(df, max_dist_between_blocks=1000000):
    temp = deepcopy(df[["start", "end", "color"]]).reset_index(drop=False)
    row_list = [list(row) for row in temp.iloc[0:1, :].itertuples(index=False)]

    for row in temp.iloc[1:, :].itertuples(index=False):
        if (row[0] == row_list[-1][0]) and (row[1] == row_list[-1][1]):
            if (row[2] - row_list[-1][3]) < max_dist_between_blocks:
                row_list[-1][3] = max(row[3], row_list[-1][3])
            else:
                row_list.append(list(row))
        else:
            row_list.append(list(row))
    # print(row_list)
    df = pd.DataFrame.from_records(row_list, columns=["scaffold", "query", "start", "end", "color"],
                                   index=["scaffold", "query"]) # index=["scaffold", "query"]
    df["target_len"] = df["end"] - df["start"]
    return df


def bed_dict_to_xlsx(bed_dict, output_prefix):
    # ---- Color configuration for xlsx file----
    # -------- Color codes ----
    green_hex = "#00FF00"
    light_green_hex = "#90EE90"
    light_blue_hex = "#ADD8E6"
    light_yellow_hex = "#F4EA56"
    light_orange_hex = "#FFD580"
    light_red_hex = "#FF7377"
    # --------
    # ----

    writer = pd.ExcelWriter('{0}.xlsx'.format(output_prefix), engine='xlsxwriter')
    workbook = writer.book
    # Adjust default format
    workbook.formats[0].set_align('center')

    long_block_format = workbook.add_format({'bg_color': light_green_hex})
    short_block_format = workbook.add_format({'bg_color': light_blue_hex})
    too_short_block_format = workbook.add_format({'bg_color': light_orange_hex})

    species_format_dict = {}
    for species in bed_dict:
        species_format_dict[species] = {}
        for scaffold in species_color_df_dict[species].index:
            species_format_dict[species][scaffold] = workbook.add_format(
                {'bg_color': species_color_df_dict[species].loc[scaffold, "color"]})

    column_start = 0

    for species in bed_dict:  # bed_col_dict:

        bed_dict[species].records["target_len"] = bed_dict[species].records["end"] - bed_dict[species].records["start"]
        if ("query_end" in bed_dict[species].records.columns) and ("query_start" in bed_dict[species].records.columns):
            bed_dict[species].records["query_len"] = bed_dict[species].records["query_end"] - bed_dict[species].records["query_start"]

        bed_dict[species].records.to_excel(writer, sheet_name=species, freeze_panes=(1, 1))
        column_number = len(bed_dict[species].records.columns) + len(bed_dict[species].records.index.names)
        row_number = len(bed_dict[species].records)

        # ----- color query and scaffold_columns -----
        scaffold_column = 0
        if "query" in bed_dict[species].records.columns:
            query_column = list(bed_dict[species].records.columns).index("query") + len(bed_dict[species].records.index.names)
        else:
            query_column = bed_dict[species].records.index.names.index("query")
        if "color" in bed_dict[species].records.columns:
            color_column = list(bed_dict[species].records.columns).index("color") + len(bed_dict[species].records.index.names)

        query_data = bed_dict[species].records["query"] if "query" in bed_dict[species].records.columns else bed_dict[species].records.index.get_level_values("query")
        for row in range(1, row_number + 1):
            writer.sheets[species].write(row, query_column, query_data[row - 1],
                                         # color query column
                                         species_format_dict[species][query_data[row - 1]])
            if "color" in bed_dict[species].records.columns:
                writer.sheets[species].write(row, color_column, bed_dict[species].records["color"][row - 1],
                                            # color color column
                                            species_format_dict[species][query_data[row - 1]])

        writer.sheets[species].set_column(column_start, len(bed_dict[species].records.columns) + len(bed_dict[species].records.index.names) - 1, 15)  #

        first_len_col = column_number - (2 if "query_len" in bed_dict[species].records.columns else 1)
        writer.sheets[species].conditional_format(1, first_len_col,
                                                  row_number, column_number - 1,
                                                  {'type': 'cell',
                                                   'criteria': 'between',
                                                   'minimum': 1000000,
                                                   'maximum': 5000000,
                                                   'format': short_block_format
                                                   })
        writer.sheets[species].conditional_format(1, first_len_col,
                                                  row_number, column_number - 1,
                                                  {'type': 'cell',
                                                   'criteria': '>=',
                                                   'value': 5000000,
                                                   'format': long_block_format
                                                   })
        writer.sheets[species].conditional_format(1, first_len_col,
                                                  row_number, column_number - 1,
                                                  {'type': 'cell',
                                                   'criteria': '<=',
                                                   'value': 1000000,
                                                   'format': too_short_block_format
                                                   })
    writer.save()


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True, type=split_comma_separated_list,
                    help="Comma-separated list of psl files with synteny blocks. "
                         "Target genome MUST be same in all files")
parser.add_argument("--input_format", action="store", dest="input_format", default="psl",
                    help="Format of input files with synteny. Allowed: psl(default), bed, bed_with_color")
parser.add_argument("--query_labels", action="store", dest="query_labels", required=True,
                    type=split_comma_separated_list,
                    help="Comma-separated list of labels for query genomes. Must follow the same order as for psl files")
parser.add_argument("--query_scaffold_white_lists", action="store", dest="query_scaffold_white_lists", default=None,
                    type=split_comma_separated_list,
                    help="Comma-separated list of files containing the only scaffolds from query genomes to include."
                         " Default: all")
parser.add_argument("--query_scaffold_black_lists", action="store", dest="query_scaffold_black_lists", default=None,
                    type=split_comma_separated_list,
                    help="Comma-separated list of files containing scaffolds from query genomes to skip at drawing. "
                         "Default: not set")

parser.add_argument("--query_color_filelist", action="store", dest="query_color_filelist", default=None,
                    type=split_comma_separated_list,
                    help="Comma-separated list of tsv files containing color sets for chromosomes of query genomes."
                         "Must follow the same order as for psl files"
                         "If not set colors will be generated automatically and differ from run to run. "
                         "If you want to draw multiple figures with the same color schemes it is recommended to run "
                         "script first time with this parameter not set and for the rest of runs use generated files"
                         " {output_prefix}.{query_label}.chr_colors.tsv. Default: not set ")
parser.add_argument("--reference_label", action="store", dest="reference_label", required=True, type=str,
                    help="Label of reference genome")
parser.add_argument("--reference_highlight_file", action="store", dest="reference_highlight_file",
                    type=lambda s: pd.read_csv(s, header=0, index_col=0, sep="\t"),
                    help="Tab-separated file with two columns ('scaffold' and 'color'). "
                         "Scaffold ids are ids after renaming"
                         "Must contain header.")
parser.add_argument("--reference_centromere_bed", action="store", dest="reference_centromere_bed", required=False,
                    type=str,
                    help="Bed file with coordinates of centromeres in reference")
parser.add_argument("--reference_scaffold_white_list", action="store", dest="reference_scaffold_white_list", default=None,
                    #type=lambda s: pd.read_csv(s, header=None, squeeze=True)) if os.path.exists(s) else pd.Series(s.split(",")),
                    type=lambda s: pd.read_csv(s, header=None).squeeze("columns") if os.path.exists(s) else pd.Series(s.split(",")),
                    help="Comma-separated list of the only scaffolds to draw (reference white list). Default: all")
parser.add_argument("--reference_scaffold_black_list", action="store", dest="reference_scaffold_black_list", default=None,
                    type=lambda s: pd.read_csv(s, header=None).squeeze("columns") if os.path.exists(s) else pd.Series(s.split(",")),
                    help="Comma-separated list of scaffolds to skip at drawing (reference black list). Default: not set")

parser.add_argument("--query_scaffold_syn_files", action="store", dest="query_scaffold_syn_files", type=split_comma_separated_list,
                    help="Comma-separated list of files with scaffold id synonyms for query genomes")
parser.add_argument("--reference_scaffold_syn_file", action="store", dest="reference_scaffold_syn_file",
                    help="File with scaffold id synonyms for reference genome")
parser.add_argument("--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in synonym file. Default: 0")
parser.add_argument("--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in synonym file synonym. Default: 1")

parser.add_argument("-z", "--reference_scaffold_order_list", action="store", dest="reference_scaffold_order_list",
                    required=True,
                    type=lambda s: pd.read_csv(s, header=None).squeeze("columns") if os.path.exists(s) else pd.Series(s.split(",")),
                    help="Comma-separated list of scaffolds to draw first and exactly in same order. "
                         "Default: not set")
parser.add_argument("--query_scaffold_order_list", action="store", dest="query_scaffold_order_lists", required=True,
                    type=split_comma_separated_list,
                    help="Comma-separated list of scaffolds to draw first and exactly in same order. "
                         "Default: not set")
parser.add_argument("-n", "--reference_scaffold_length_file", action="store", dest="reference_scaffold_length_file",
                    required=True,
                    help="File with lengths of scaffolds for reference genome")
parser.add_argument("--initial_min_block_len_list", action="store", dest="initial_min_block_len_list",
                    type=lambda s: list(map(int, s.split(","))), default=[0,],
                    help="Comma-separated list of minimal block length for initial filtration. Default: 0")

parser.add_argument("--secondary_min_block_len_list", action="store", dest="secondary_min_block_len_list",
                    type=lambda s: list(map(int, s.split(","))), default=[1000000, ],
                    help="Comma-separated list of minimal block length for secondary filtration. Default: 1000000")

parser.add_argument("--max_dist_between_short_blocks_list", action="store", dest="max_dist_between_short_blocks_list",
                    type=lambda s: list(map(int, s.split(","))), default=[3000000, ],
                    help="Comma-separated list of maximal distance between short blocks "
                         "during second stage of filtration"
                         "Blocks being shorter than --secondary_min_block_len_list and having both upstream and "
                         "downstream distance for same chromosome (both target and query) blocks will be discarded."
                         "Default: 3000000")
parser.add_argument("--max_dist_between_blocks_list", action="store", dest="max_dist_between_blocks_list",
                    type=lambda s: sorted(list(map(int, s.split(",")))), default=[1000000, ],
                    help="Comma-separated list of maximal distance between blocks for merging of adjacent blocks"
                         "Default: 1000000")
parser.add_argument("--final_min_block_len_list", action="store", dest="final_min_block_len_list",
                    type=lambda s: sorted(list(map(int, s.split(",")))), default=[1000000, ],
                    help="Comma-separated list of minimal block length for final filtration."
                         "Default: 1000000")


parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=split_comma_separated_list,
                    default=("png", "svg"),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg,png")

parser.add_argument("-l", "--title", action="store", dest="title", default="Coverage",
                    help="Suptitle of figure. Default: 'Synteny'")
parser.add_argument("--title_fontsize", action="store", dest="title_fontsize", default=20, type=int,
                    help="Fontsize of the figure. Default: 20")

parser.add_argument("--stranded", action="store_true", dest="stranded", default=False,
                    help="Stranded features and tracks. Default: False")
parser.add_argument("--rounded", action="store_true", dest="rounded", default=False,
                    help="Rounded tracks. Default: False")
parser.add_argument("--stranded_end", action="store_true", dest="stranded_end", default=False,
                    help="Stranded ends for tracks. Works only if --stranded is set. Default: False")
parser.add_argument("--hide_track_label", action="store_true", dest="hide_track_label", default=False,
                    help="Hide track label. Default: False")
parser.add_argument("--feature_shape", action="store", dest="feature_shape", default="rectangle",
                    help="Shape of features. Allowed: rectangle(default), circle, ellipse")

parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float, default=0.2,
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
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Print additional info to stdout")

parser.add_argument("--subplot_scale", action="store_true", dest="subplot_scale",
                    help="Scale feature x size by subplot x/y ratio. Default: off")
parser.add_argument("--track_group_scale", action="store_true", dest="track_group_scale",
                    help="Scale feature x size by track_group x/y ratio. Default: off")

parser.add_argument("--x_tick_fontsize", action="store", dest="x_tick_fontsize", type=int, default=None,
                    help="Fontsize of xticks. Default: matplotlib default")
parser.add_argument("--invert_coordinates_for_target_negative_strand", action="store_true",
                    dest="invert_coordinates_for_target_negative_strand",
                    default=False,
                    help="Invert coordinates for target negative strand. Default: False")


args = parser.parse_args()


reference = args.reference_label
query_list = args.query_labels

species_chr_syn_dict = {query: pd.read_csv(syn_file,
                                           usecols=(args.syn_file_key_column, args.syn_file_value_column),
                                           header=None,
                                           sep="\t", index_col=args.syn_file_key_column).squeeze("columns") if syn_file else None for query, syn_file in zip(query_list, args.query_scaffold_syn_files)}

species_chr_syn_dict[reference] = pd.read_csv(args.reference_scaffold_syn_file,
                                              usecols=(args.syn_file_key_column,args.syn_file_value_column,),
                                              header=None,
                                              sep="\t", index_col=args.syn_file_key_column).squeeze("columns") if args.reference_scaffold_syn_file else None

if args.reference_centromere_bed:
    centromere_df = pd.read_csv(args.reference_centromere_bed,
                                usecols=(0, 1, 2),
                                index_col=0,
                                header=None,
                                sep="\t", names=["scaffold_id", "start", "end"])
    centromere_df.rename(index=species_chr_syn_dict[reference], inplace=True)
    #print(centromere_df)
else:
    centromere_df = None

if args.query_scaffold_white_lists:
    species_white_list_dict = {query: pd.read_csv(white_list_file, header=None).squeeze("columns") if white_list_file else None for query, white_list_file in zip(query_list, args.query_scaffold_white_lists)}
else:
    species_white_list_dict = {query: None for query in query_list}

species_white_list_dict[reference] = args.reference_scaffold_white_list #pd.read_csv(args.reference_scaffold_white_list, header=None, squeeze=True) if args.reference_scaffold_white_list else None

if args.query_scaffold_black_lists:
    species_black_list_dict = {query: pd.read_csv(black_list_file, header=None).squeeze("columns") if black_list_file else None for query, black_list_file in zip(query_list, args.query_scaffold_black_lists)}
else:
    species_black_list_dict = {query: None for query in query_list}

species_black_list_dict[reference] = args.reference_scaffold_black_list #pd.read_csv(args.reference_scaffold_black_list, header=None, squeeze=True) if args.reference_scaffold_black_list else None
species_orderlist_dict = {query: pd.read_csv(orderlist_file, header=None).squeeze("columns").iloc[::-1] for query, orderlist_file in zip (query_list, args.query_scaffold_order_lists)}
species_orderlist_dict[reference] = args.reference_scaffold_order_list[::-1]

species_color_df_dict = OrderedDict()
if args.query_color_filelist is None:
    color_number = max([len(species_orderlist_dict[query]) for query in query_list])
    colors = distinctipy.get_colors(color_number)
    color_list = list(map(rgb_tuple_to_hex, colors))

    for species in query_list:
        species_color_df_dict[species] = pd.DataFrame()
        species_color_df_dict[species]["scaffold"] = species_orderlist_dict[species]
        species_color_df_dict[species]["color"] = color_list[:len(species_orderlist_dict[species])]
        species_color_df_dict[species].set_index("scaffold", inplace=True)
        species_color_df_dict[species].to_csv("{}.{}.chr_colors.tsv".format(args.output_prefix, species), sep="\t", header=True, index=True)

else:
    for species, color_file in zip(query_list, args.query_color_filelist):
        species_color_df_dict[species] = pd.read_csv(color_file, sep="\t", index_col=0, header=0)
        species_color_df_dict[species].to_csv("{}.{}.chr_colors.tsv".format(args.output_prefix, species), sep="\t",
                                              header=True, index=True)

reference_scaffold_length_df = pd.read_csv(args.reference_scaffold_length_file, sep='\t',
                                           header=None, names=("scaffold", "length"), index_col=0)
reference_scaffold_length_df.index = pd.Index(list(map(str, reference_scaffold_length_df.index)))
if not species_chr_syn_dict[reference].empty:
    reference_scaffold_length_df.rename(index=species_chr_syn_dict[reference], inplace=True)

bed_col_dict = OrderedDict()

if args.input_format == "psl":
    psl_col_dict = {query: CollectionPSL(in_file=psl, parsing_mode="coordinates_only",
                                         target_syn_dict=species_chr_syn_dict[reference],
                                         target_black_list=species_black_list_dict[reference],
                                         target_white_list=species_white_list_dict[reference],
                                         query_syn_dict=species_chr_syn_dict[query],
                                         query_black_list=species_black_list_dict[query],
                                         query_white_list=species_white_list_dict[query],
                                         invert_coordinates_for_target_negative_strand=args.invert_coordinates_for_target_negative_strand
                                         ) for query, psl in zip(query_list, args.input)}

    for species in query_list:
        bed_col_dict[species] = CollectionBED(
            records=psl_col_dict[species].records[["tName", "tStart", "tEnd", "qName", "qStart", "qEnd", "strand"]],
            records_columns=["scaffold", "start", "end", "query", "query_start", "query_end", "strand"],
            format="bed_synteny_track",
            parsing_mode="all")
        bed_col_dict[species].records.set_index("scaffold", inplace=True)

        bed_col_dict[species].records["color"] = bed_col_dict[species].records["query"].replace(species_color_df_dict[species]["color"])

        bed_col_dict[species].records.sort_values(by=["scaffold", "start", "end", ], inplace=True)

elif args.input_format in ["bed", "bed_with_color"]:
    bed_file_dict = {query: bed for query, bed in zip(query_list, args.input)}

    for species in query_list:
        bed_col_dict[species] = CollectionBED(in_file=bed_file_dict[species], header_in_file=True,
                                              format="bed_synteny_track", parsing_mode="all",
                                              scaffold_syn_dict=species_chr_syn_dict[reference] if species_chr_syn_dict[
                                                                                                       reference] is not None else None,
                                              rename_dict={"query": species_chr_syn_dict[species]} if
                                              species_chr_syn_dict[species] is not None else None)

        bed_col_dict[species].records.sort_values(by=["scaffold", "start", "end", ], inplace=True)
        if args.input_format != "bed_with_color":
            bed_col_dict[species].records["color"] = bed_col_dict[species].records["query"].replace(species_color_df_dict[species]["color"])
else:
    raise ValueError("ERROR!!! Unrecognized format of the input file(s)!")
query_species_color_df_dict = {sp: species_color_df_dict[sp] for sp in query_list}

#for query in bed_col_dict:
#    print(bed_col_dict[query].records)
#print(bed_col_dict[query_list[0]].records)


for min_block_length in args.initial_min_block_len_list:
    prefiltered_bed_col_dict = {}
    if min_block_length == 0:
        prefiltered_bed_col_dict = bed_col_dict
    else:
        prefiltered_bed_col_dict = deepcopy(bed_col_dict)
        for species in prefiltered_bed_col_dict:
            prefiltered_bed_col_dict[species].records = prefiltered_bed_col_dict[species].records[(prefiltered_bed_col_dict[species].records["end"] - prefiltered_bed_col_dict[species].records["start"]) >= min_block_length]

    for species in prefiltered_bed_col_dict:  # bed_col_dict:
        # ---- save original blocks to bed ----
        prefiltered_bed_col_dict[species].records.to_csv("{0}.{1}.to.{2}.initial_min_block_len_{3}.tsv".format(args.output_prefix,
                                                                                                    species,
                                                                                                    reference,
                                                                                                    min_block_length),
                                                      sep="\t",
                                                      header=True, index=True)
        # -----

    Visualization.draw_features(prefiltered_bed_col_dict,
                                reference_scaffold_length_df,
                                species_orderlist_dict[reference],
                                "{0}.initial_min_block_len_{1}".format(args.output_prefix, min_block_length),
                                legend=Visualization.chromosome_legend(query_species_color_df_dict,
                                                                       species_orderlist_dict[reference]),
                                centromere_df=centromere_df,
                                highlight_df=args.reference_highlight_file,
                                figure_width=args.figure_width, figure_height_per_scaffold=args.figure_height_per_scaffold,
                                dpi=300,
                                #colormap=None, thresholds=None, colors=None, background=None,
                                default_color="red",  # TODO: check if it is possible to remove it
                                title=args.title,
                                extensions=args.output_formats,
                                feature_start_column_id="start",
                                feature_end_column_id="end",
                                feature_color_column_id="color",
                                feature_length_column_id="length",
                                feature_height_fraction=0.7,
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
    bed_dict_to_xlsx(prefiltered_bed_col_dict, '{0}.initial_min_block_len_{1}'.format(args.output_prefix, min_block_length))

    """
    writer = pd.ExcelWriter('{0}.initial_min_block_len_{1}.xlsx'.format(args.output_prefix, min_block_length),
                            engine='xlsxwriter')
    workbook = writer.book
    # Adjust default format
    workbook.formats[0].set_align('center')

    long_block_format = workbook.add_format({'bg_color': light_green_hex})
    short_block_format = workbook.add_format({'bg_color': light_blue_hex})
    too_short_block_format = workbook.add_format({'bg_color': light_orange_hex})

    species_format_dict = {}
    for species in prefiltered_bed_col_dict:
        species_format_dict[species] = {}
        for scaffold in species_color_df_dict[species].index:
            species_format_dict[species][scaffold] = workbook.add_format({'bg_color': species_color_df_dict[species].loc[scaffold, "color"]})

    column_start = 0

    for species in prefiltered_bed_col_dict:  # bed_col_dict:

        prefiltered_bed_col_dict[species].records["target_len"] = prefiltered_bed_col_dict[species].records["end"] - prefiltered_bed_col_dict[species].records["start"]
        prefiltered_bed_col_dict[species].records["query_len"] = prefiltered_bed_col_dict[species].records["query_end"] - prefiltered_bed_col_dict[species].records["query_start"]

        prefiltered_bed_col_dict[species].records.to_excel(writer, sheet_name=species, freeze_panes=(1, 1))
        column_number = len(prefiltered_bed_col_dict[species].records.columns)
        row_number = len(prefiltered_bed_col_dict[species].records)

        # ----- color query and scaffold_columns -----
        scaffold_column = 0
        query_column = list(prefiltered_bed_col_dict[species].records.columns).index("query") + 1  # +1 added to take into account one-level index
        color_column = list(prefiltered_bed_col_dict[species].records.columns).index("color") + 1  # +1 added to take into account one-level index

        for row in range(1, row_number + 1):
            writer.sheets[species].write(row, query_column, prefiltered_bed_col_dict[species].records["query"][row-1], # color query column
                                         species_format_dict[species][prefiltered_bed_col_dict[species].records["query"][row-1]])
            writer.sheets[species].write(row, color_column, prefiltered_bed_col_dict[species].records["color"][row - 1], # color color column
                                         species_format_dict[species][prefiltered_bed_col_dict[species].records["query"][row - 1]])

        writer.sheets[species].set_column(column_start, len(prefiltered_bed_col_dict[species].records.columns), 15)  #

        writer.sheets[species].conditional_format(1, column_number - 2,
                                                  row_number, column_number,
                                                  {'type': 'cell',
                                                   'criteria': 'between',
                                                   'minimum': 1000000,
                                                   'maximum': 5000000,
                                                   'format': short_block_format
                                                   })
        writer.sheets[species].conditional_format(1, column_number - 2,
                                                  row_number, column_number,
                                                  {'type': 'cell',
                                                   'criteria': '>=',
                                                   'value': 5000000,
                                                   'format': long_block_format
                                                   })
        writer.sheets[species].conditional_format(1, column_number - 2,
                                                  row_number, column_number,
                                                  {'type': 'cell',
                                                   'criteria': '<=',
                                                   'value': 1000000,
                                                   'format': too_short_block_format
                                                   })
    writer.save()
    """

    for secondary_min_block_len in args.secondary_min_block_len_list:
        for max_dist_between_short_blocks in args.max_dist_between_short_blocks_list:
            filtered_bed_col_dict = deepcopy(prefiltered_bed_col_dict)

            for species in filtered_bed_col_dict:
                filtered_bed_col_dict[species].records = filter_isolated_short_blocks(filtered_bed_col_dict[species].records,
                                                                                      min_block_len=secondary_min_block_len,
                                                                                      max_dist_between_short_blocks=max_dist_between_short_blocks)

            second_stage_output_suffix = "initial_min_block_len_{0}.secondary_min_block_len_{1}.max_dist_between_short_blocks_{2}".format(min_block_length,
                                                                                                                                          secondary_min_block_len,
                                                                                                                                          max_dist_between_short_blocks)

            for species in filtered_bed_col_dict:  # bed_col_dict:
                # ---- save original blocks to bed ----
                filtered_bed_col_dict[species].records.to_csv(
                    "{0}.{1}.to.{2}.{3}.tsv".format(args.output_prefix,
                                                    species,
                                                    reference,
                                                    second_stage_output_suffix),
                    sep="\t",
                    header=True, index=True)
                # -----

                Visualization.draw_features(filtered_bed_col_dict,
                                            reference_scaffold_length_df,
                                            species_orderlist_dict[reference],
                                            "{0}.{1}".format(args.output_prefix,second_stage_output_suffix),
                                            legend=Visualization.chromosome_legend(query_species_color_df_dict,
                                                                                   species_orderlist_dict[reference]),
                                            centromere_df=centromere_df,
                                            highlight_df=args.reference_highlight_file,
                                            figure_width=args.figure_width,
                                            figure_height_per_scaffold=args.figure_height_per_scaffold,
                                            dpi=300,
                                            # colormap=None, thresholds=None, colors=None, background=None,
                                            default_color="red",  # TODO: check if it is possible to remove it
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
                bed_dict_to_xlsx(filtered_bed_col_dict,
                                 "{0}.{1}".format(args.output_prefix,
                                                  second_stage_output_suffix))

            for max_dist_between_blocks in args.max_dist_between_blocks_list:
                for species in filtered_bed_col_dict:
                    filtered_bed_col_dict[species].records = merge_adjacent_blocks(filtered_bed_col_dict[species].records,
                                                                                   max_dist_between_blocks=max_dist_between_blocks)
                    #filtered_bed_col_dict[species].records["color"] =
                third_stage_output_suffix = "max_dist_between_adjacent_blocks_{0}".format(max_dist_between_blocks)

                for species in filtered_bed_col_dict:  # bed_col_dict:
                    # ---- save original blocks to bed ----
                    filtered_bed_col_dict[species].records.to_csv("{0}.{1}.to.{2}.{3}.{4}.tsv".format(args.output_prefix,
                                                                                                      species,
                                                                                                      reference,
                                                                                                      second_stage_output_suffix,
                                                                                                      third_stage_output_suffix),
                                                                  sep="\t",
                                                                  header=True, index=True)
                    # -----

                    Visualization.draw_features(filtered_bed_col_dict,
                                                reference_scaffold_length_df,
                                                species_orderlist_dict[reference],
                                                "{0}.{1}.{2}".format(args.output_prefix,
                                                                     second_stage_output_suffix,
                                                                     third_stage_output_suffix),
                                                legend=Visualization.chromosome_legend(query_species_color_df_dict,
                                                                                       species_orderlist_dict[
                                                                                           reference]),
                                                centromere_df=centromere_df,
                                                highlight_df=args.reference_highlight_file,
                                                figure_width=args.figure_width,
                                                figure_height_per_scaffold=args.figure_height_per_scaffold,
                                                dpi=300,
                                                # colormap=None, thresholds=None, colors=None, background=None,
                                                default_color="red",  # TODO: check if it is possible to remove it
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
                                                stranded_tracks=False,
                                                rounded_tracks=args.rounded,
                                                stranded_end_tracks=False,
                                                xtick_fontsize=args.x_tick_fontsize,
                                                subplot_title_fontsize=args.title_fontsize,
                                                subplot_title_fontweight='bold'
                                                )
                    bed_dict_to_xlsx(filtered_bed_col_dict,
                                     "{0}.{1}.{2}".format(args.output_prefix,
                                                          second_stage_output_suffix,
                                                          third_stage_output_suffix))

                for final_min_block_len in args.final_min_block_len_list:
                    for species in filtered_bed_col_dict:
                        filtered_bed_col_dict[species].records = filtered_bed_col_dict[species].records[filtered_bed_col_dict[species].records["end"] - filtered_bed_col_dict[species].records["start"] > final_min_block_len]

                    forth_stage_output_suffix = "final_min_block_len_{0}".format(final_min_block_len)

                    for species in filtered_bed_col_dict:  # bed_col_dict:
                        # ---- save original blocks to bed ----
                        filtered_bed_col_dict[species].records.to_csv(
                            "{0}.{1}.to.{2}.{3}.{4}.{5}.tsv".format(args.output_prefix,
                                                                    species,
                                                                    reference,
                                                                    second_stage_output_suffix,
                                                                    third_stage_output_suffix,
                                                                    forth_stage_output_suffix),
                            sep="\t",
                            header=True, index=True)
                        # -----

                        Visualization.draw_features(filtered_bed_col_dict,
                                                    reference_scaffold_length_df,
                                                    species_orderlist_dict[reference],
                                                    "{0}.{1}.{2}.{3}".format(args.output_prefix,
                                                                             second_stage_output_suffix,
                                                                             third_stage_output_suffix,
                                                                             forth_stage_output_suffix),
                                                    legend=Visualization.chromosome_legend(query_species_color_df_dict,
                                                                                           species_orderlist_dict[
                                                                                               reference]),
                                                    centromere_df=centromere_df,
                                                    highlight_df=args.reference_highlight_file,
                                                    figure_width=args.figure_width,
                                                    figure_height_per_scaffold=args.figure_height_per_scaffold,
                                                    dpi=300,
                                                    # colormap=None, thresholds=None, colors=None, background=None,
                                                    default_color="red",  # TODO: check if it is possible to remove it
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
                                                    stranded_tracks=False,
                                                    rounded_tracks=args.rounded,
                                                    stranded_end_tracks=False,
                                                    xtick_fontsize=args.x_tick_fontsize,
                                                    subplot_title_fontsize=args.title_fontsize,
                                                    subplot_title_fontweight='bold'
                                                    )
                        bed_dict_to_xlsx(filtered_bed_col_dict,
                                         "{0}.{1}.{2}.{3}".format(args.output_prefix,
                                                                  second_stage_output_suffix,
                                                                  third_stage_output_suffix,
                                                                  forth_stage_output_suffix))
