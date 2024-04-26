#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import glob
import argparse
from copy import deepcopy
from pathlib import Path
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


def expand_path(path_template: str, skip=False):
    if (path_template[0] == "/") or skip:
        print("Skipping expanding path {0} as it is global path or expanding was turned off ...".format(path_template))
        # avoid expanding global paths
        return path_template

    path_list = list(Path("./").glob(path_template))
    if len(path_list) > 1:
        raise ValueError("ERROR!!! There is more than one file corresponding to the template {0} ...".format(path_template))
    elif len(path_list) == 0:
        raise ValueError("ERROR!!! There are no files corresponding to the template {0} ...".format(path_template))
    else:
        return str(path_list[0])


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
    #print(tmp_df.groupby(["scaffold",
    #                         "query"], sort=False, group_keys=False).apply(dist_to_next_seg_from_same_chr))
    tmp_df = tmp_df.groupby(["scaffold",
                             "query"], sort=False, group_keys=False).apply(dist_to_next_seg_from_same_chr).sort_values(by=["scaffold",
                                                                                                         "start",
                                                                                                         "end",
                                                                                                         "query",
                                                                                                         "query_start",
                                                                                                         "query_end"],
                                                                                                     axis=0)  # .reset_index(level=(0,1), drop=True)
    tmp_df["embedded"].fillna(False, inplace=True)
    # recalculate distances after removal of embedded blocks
    tmp_df = tmp_df[tmp_df["embedded"] <= 0]
    tmp_df = tmp_df.groupby(["scaffold",
                             "query"], sort=False, group_keys=False).apply(dist_to_next_seg_from_same_chr).sort_values(by=["scaffold",
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
    #print(bed_dict)
    for species in bed_dict:
        species_format_dict[species] = {}
        for scaffold in color_df_dict[species].index:
            species_format_dict[species][scaffold] = workbook.add_format(
                {'bg_color': color_df_dict[species].loc[scaffold, "color"]})

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

        query_data = list(bed_dict[species].records["query"] if "query" in bed_dict[species].records.columns else bed_dict[species].records.index.get_level_values("query"))
        #print(query_data)
        for row in range(1, row_number + 1):
            writer.sheets[species].write(row, query_column, query_data[row - 1],
                                         # color query column
                                         species_format_dict[species][query_data[row - 1]])
            if "color" in bed_dict[species].records.columns:
                writer.sheets[species].write(row, color_column, bed_dict[species].records["color"].iloc[row - 1],
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
    writer.close()


def get_filenames_for_extension(dir_path, extension_list, force_uniq=True):
    filelist = []
    for extension in extension_list:
        filelist += list(glob.glob(str(dir_path) + "/*{0}".format(extension)))
        #print(extension, filelist)
    if not filelist:
        return None
    print(filelist)
    if force_uniq:
        if len(filelist) > 1:
            raise ValueError("Found more than one file with extensions: {0} in directory {1}".format(",".join(extension_list, str(dir_path))))
        else:
            return filelist[0]

    return filelist


def invert_coordinates_in_synteny_table(df, scaffold_list, length_df, scaffold_column, start_column, end_column, strand_column, inverted_scaffolds_label="'"):
    temp_df = deepcopy(df)
    temp_df["length_column"] = 0
    original_index_name = temp_df.index.name
    columns_list = list(temp_df.columns)
    temp_df.reset_index(drop=False, inplace=True)
    temp_df.set_index(scaffold_column, inplace=True)
    #temp_df.to_csv("tmp", sep="\t", index=True, header=True)
    for scaffold in temp_df.index.unique():
        temp_df.loc[scaffold, "length_column"] = length_df.loc[scaffold, "length"]

    temp_df.loc[temp_df.index.isin(scaffold_list), start_column], temp_df.loc[temp_df.index.isin(scaffold_list), end_column] = temp_df.loc[temp_df.index.isin(scaffold_list), "length_column"] - temp_df.loc[temp_df.index.isin(scaffold_list), end_column], \
               temp_df.loc[temp_df.index.isin(scaffold_list), "length_column"] - temp_df.loc[temp_df.index.isin(scaffold_list), start_column]

    plus_indexes, minus_indexes = (temp_df[strand_column] == "+") & temp_df.index.isin(scaffold_list), (temp_df[strand_column] == "-") & temp_df.index.isin(scaffold_list)
    temp_df.loc[plus_indexes, strand_column], temp_df.loc[minus_indexes, strand_column] = "-", "+"
    temp_df.reset_index(drop=False, inplace=True)
    if inverted_scaffolds_label is not None:
        for scaffold in scaffold_list:
            temp_df.loc[temp_df[scaffold_column] == scaffold, scaffold_column] = scaffold + inverted_scaffolds_label
    if (original_index_name is not None) and original_index_name in temp_df.columns:
        temp_df.set_index(original_index_name, inplace=True)

    return temp_df[columns_list]  # remove added length column and restore column order


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True,
                    help="Input directory with data. Must contain subfolder for each genome. "
                         "Subfolders should have same name as genomes in --genome_orderlist"
                         "Each subfolder should contain: *.whitelist, *.len and synteny file "
                         "(except for the last genome). *.orderlist, *.invertlist and *.syn file are optional")
parser.add_argument("--synteny_format", action="store", dest="synteny_format", default="psl",
                    help="Format of the synteny file. Allowed: psl(default), bed, bed_with_color")
parser.add_argument("--query_orderlist", action="store", dest="query_orderlist", required=True,
                    type=split_comma_separated_list,
                    help="Comma-separated list of labels for query genomes.")
parser.add_argument("--invert_genome_order", action="store_true", dest="invert_genome_order", default=False,
                    help="Invert order of the genomes in the --genome_orderlist. Default: False")

parser.add_argument("--use_original_colors", action="store_true", dest="use_original_colors", default=False,
                    help="Use colors from .color file. If not set colors will be select ")
parser.add_argument("--reference_label", action="store", dest="reference_label", required=True, type=str,
                    help="Label of reference genome")
parser.add_argument("--reference_highlight_file", action="store", dest="reference_highlight_file",
                    type=lambda s: pd.read_csv(s, header=0, index_col=0, sep="\t"),
                    help="Tab-separated file with two columns ('scaffold' and 'color'). "
                         "Scaffold ids are ids after renaming"
                         "Must contain header.")
#parser.add_argument("--reference_centromere_bed", action="store", dest="reference_centromere_bed", required=False,
#                    type=str,
#                    help="Bed file with coordinates of centromeres in reference")
parser.add_argument("--reference_scaffold_white_list", action="store", dest="reference_scaffold_white_list", default=None,
                    #type=lambda s: pd.read_csv(s, header=None, squeeze=True)) if os.path.exists(s) else pd.Series(s.split(",")),
                    type=lambda s: pd.read_csv(s, header=None).squeeze("columns") if os.path.exists(s) else pd.Series(s.split(",")),
                    help="Comma-separated list of the only scaffolds to draw (reference white list). Default: use  .whitelist from reference genome folder")
parser.add_argument("-z", "--reference_scaffold_order_list", action="store", dest="reference_scaffold_order_list",
                    type=lambda s: pd.read_csv(s, header=None).squeeze("columns") if os.path.exists(s) else pd.Series(s.split(",")),
                    help="Comma-separated list of scaffolds to draw first and exactly in same order. "
                         "Default: not set")
parser.add_argument("--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in synonym file. Default: 0")
parser.add_argument("--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in synonym file synonym. Default: 1")
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

parser.add_argument("-l", "--title", action="store", dest="title", default="Synteny",
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
parser.add_argument("--ymax_multiplier", action="store", dest="ymax_multiplier", type=float, default=1.0,
                    help="Multiplier for y max limit of figure. Default: 1.0")
parser.add_argument("--figure_header_height", action="store", dest="figure_header_height", type=int, default=0,
                    help="Additional height of figure to account for header. Default: 0")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Print additional info to stdout")
parser.add_argument("--hide_legend", action="store_true", dest="hide_legend", default=False,
                    help="Don't draw legend. Default: False")
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
parser.add_argument("-x", "--expand_paths", action="store_true", dest="expand_paths", default=False,
                    help="Expand paths by following regular expressions")

args = parser.parse_args()


reference = args.reference_label
query_list = args.query_orderlist

data_dir = args.input_dir
data_dir_path = Path(data_dir)

genome_list = query_list + [reference]
syn_file_key_column, syn_file_value_column = args.syn_file_key_column, args.syn_file_value_column

synteny_format = args.synteny_format

# read files
#print(data_dir_path)
#print(genome_list)
# whitelist is obligatory
for genome in genome_list:
    print(data_dir_path / genome)
    #print(get_filenames_for_extension(data_dir_path / genome, extension_list=["whitelist"]))
#print(data_dir_path / genome)
#print(get_filenames_for_extension(data_dir_path / genome, extension_list=["whitelist"]))


whitelist_series_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["whitelist"]),
                                             sep="\t", header=None, comment="#").squeeze("columns") for genome in genome_list}
if args.reference_scaffold_white_list is not None:
    whitelist_series_dict[reference] = args.reference_scaffold_white_list

# orderlist might be absent in folders
orderlist_series_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["orderlist"]),
                                             sep="\t", header=None, comment="#").squeeze("columns") if get_filenames_for_extension(data_dir_path / genome, extension_list=["orderlist"]) is not None else pd.Series(dtype=str) for genome in genome_list}
if args.reference_scaffold_order_list is not None:
    orderlist_series_dict[reference] = args.reference_scaffold_order_list
orderlist_series_dict[reference] = orderlist_series_dict[reference][::-1] if not args.invert_genome_order else orderlist_series_dict[reference]
# invertlist might be absent in folders
invertlist_series_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["invertlist"]),
                                             sep="\t", header=None, comment="#").squeeze("columns") if get_filenames_for_extension(data_dir_path / genome, extension_list=["invertlist"]) is not None else pd.Series(dtype=str) for genome in genome_list}
# lenlist is obligatory
lenlist_df_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["len"]),
                                             sep="\t", header=None, comment="#", names=["scaffold", "length"], index_col=0) for genome in genome_list}
# synfile might be absent
syn_df_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["syn"]), usecols=(syn_file_key_column, syn_file_value_column),
                                             sep="\t", header=None, comment="#",
                                   names=["key", "syn"] if syn_file_key_column <= syn_file_value_column else ["syn", "key"]).set_index("key") if get_filenames_for_extension(data_dir_path / genome, extension_list=["syn"]) is not None else pd.DataFrame(columns=["syn"]) for genome in genome_list}

color_df_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["color"]),
                                             sep="\t", header=0, names=["scaffold", "color"], index_col=0) if get_filenames_for_extension(data_dir_path / genome, extension_list=["color"]) is not None else pd.Series(dtype=str) for genome in genome_list}


for genome in genome_list:
    if not syn_df_dict[genome].empty:
        lenlist_df_dict[genome].rename(index=syn_df_dict[genome]["syn"].to_dict(), inplace=True)


"""
species_chr_syn_dict = {query: pd.read_csv(expand_path(syn_file, skip=not args.expand_paths),
                                           usecols=(args.syn_file_key_column, args.syn_file_value_column),
                                           header=None,
                                           sep="\t", index_col=args.syn_file_key_column).squeeze("columns") if syn_file else None for query, syn_file in zip(query_list, args.query_scaffold_syn_files)}

species_chr_syn_dict[reference] = pd.read_csv(expand_path(args.reference_scaffold_syn_file, skip=not args.expand_paths),
                                              usecols=(args.syn_file_key_column, args.syn_file_value_column,),
                                              header=None,
                                              sep="\t", index_col=args.syn_file_key_column).squeeze("columns") if args.reference_scaffold_syn_file else None
"""

if get_filenames_for_extension(data_dir_path / reference, extension_list=["centromere.bed"]) is None:
    centromere_df = None
else:
    centromere_df = pd.read_csv(get_filenames_for_extension(data_dir_path / reference, extension_list=["centromere.bed"]),
                                usecols=(0, 1, 2),
                                index_col=0,
                                header=None,
                                sep="\t", names=["scaffold_id", "start", "end"])
    centromere_df.rename(index=syn_df_dict[reference]["syn"].to_dict(), inplace=True)

"""
if args.reference_centromere_bed:
    centromere_df = pd.read_csv(expand_path(args.reference_centromere_bed, skip=not args.expand_paths),
                                usecols=(0, 1, 2),
                                index_col=0,
                                header=None,
                                sep="\t", names=["scaffold_id", "start", "end"])
    centromere_df.rename(index=syn_df_dict[reference]["syn"].to_dict(), inplace=True)

else:
    centromere_df = None
"""
print("Coordinates of centromere:")
print(centromere_df)
"""
if args.query_scaffold_white_lists:
    species_white_list_dict = {query: pd.read_csv(expand_path(white_list_file, skip=not args.expand_paths),
                                                  header=None).squeeze("columns") if white_list_file else None for query, white_list_file in zip(query_list, args.query_scaffold_white_lists)}
else:
    species_white_list_dict = {query: None for query in query_list}

species_white_list_dict[reference] = args.reference_scaffold_white_list #pd.read_csv(args.reference_scaffold_white_list, header=None, squeeze=True) if args.reference_scaffold_white_list else None

if args.query_scaffold_black_lists:
    species_black_list_dict = {query: pd.read_csv(expand_path(black_list_file, skip=not args.expand_paths),
                                                  header=None).squeeze("columns") if black_list_file else None for query, black_list_file in zip(query_list, args.query_scaffold_black_lists)}
else:
    species_black_list_dict = {query: None for query in query_list}

species_black_list_dict[reference] = args.reference_scaffold_black_list #pd.read_csv(args.reference_scaffold_black_list, header=None, squeeze=True) if args.reference_scaffold_black_list else None
species_orderlist_dict = {query: pd.read_csv(expand_path(orderlist_file, skip=not args.expand_paths),
                                             header=None).squeeze("columns").iloc[::-1] for query, orderlist_file in zip (query_list, args.query_scaffold_order_lists)}
species_orderlist_dict[reference] = args.reference_scaffold_order_list[::-1]
"""
color_df_dict = OrderedDict()

if not args.use_original_colors:
    color_number = max([len(orderlist_series_dict[query]) for query in query_list])
    colors = distinctipy.get_colors(color_number)
    color_list = list(map(rgb_tuple_to_hex, colors))

    for species in query_list:
        color_df_dict[species] = pd.DataFrame()
        color_df_dict[species]["scaffold"] = orderlist_series_dict[species]
        color_df_dict[species]["color"] = color_list[:len(orderlist_series_dict[species])]
        color_df_dict[species].set_index("scaffold", inplace=True)
        color_df_dict[species].to_csv("{}.{}.chr_colors.tsv".format(args.output_prefix, species), sep="\t", header=True, index=True)

else:
    color_df_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["colors"]),
                                             sep="\t", header=0, names=["scaffold", "color"], index_col=0) for genome in query_list}
    #for species, color_file in zip(query_list, args.query_color_filelist):
    #    #print(pd.read_csv(expand_path(color_file, skip=not args.expand_paths), sep="\t", index_col=0, header=0))
    #    print(species)
    #    species_color_df_dict[species] = pd.read_csv(expand_path(color_file, skip=not args.expand_paths), sep="\t", index_col=0, header=0)
    #    species_color_df_dict[species].to_csv("{}.{}.chr_colors.tsv".format(args.output_prefix, species), sep="\t",
    #                                          header=True, index=True)

#reference_scaffold_length_df = pd.read_csv(expand_path(args.reference_scaffold_length_file, skip=not args.expand_paths), sep='\t',
#                                           header=None, names=("scaffold", "length"), index_col=0)
#reference_scaffold_length_df.index = pd.Index(list(map(str, reference_scaffold_length_df.index)))
#if not species_chr_syn_dict[reference].empty:
#    reference_scaffold_length_df.rename(index=species_chr_syn_dict[reference], inplace=True)

bed_col_dict = OrderedDict()

if args.synteny_format == "psl":
    for query in query_list:
        print(get_filenames_for_extension(data_dir_path / query, extension_list=["psl", "psl.gz"]),)
    psl_col_dict = {query: CollectionPSL(in_file=get_filenames_for_extension(data_dir_path / query,
                                                                             extension_list=["psl", "psl.gz"]),
                                         parsing_mode="coordinates_only",
                                         #target_syn_dict=syn_df_dict[reference].to_dict(),
                                         target_black_list=None,
                                         target_white_list=whitelist_series_dict[reference],
                                         #query_syn_dict=syn_df_dict[query].to_dict(),
                                         query_black_list=None,
                                         query_white_list=whitelist_series_dict[query],
                                         invert_coordinates_for_target_negative_strand=args.invert_coordinates_for_target_negative_strand
                                         ) for query in query_list}
    for query in query_list:
        #print(syn_df_dict[reference])
        psl_col_dict[query].records["tName"].replace(syn_df_dict[reference]["syn"].to_dict(), inplace=True)
        psl_col_dict[query].records["qName"].replace(syn_df_dict[query]["syn"].to_dict(), inplace=True)

    for species in query_list:
        bed_col_dict[species] = CollectionBED(
            records=psl_col_dict[species].records[["tName", "tStart", "tEnd", "qName", "qStart", "qEnd", "strand"]],
            records_columns=["scaffold", "start", "end", "query", "query_start", "query_end", "strand"],
            format="bed_synteny_track",
            parsing_mode="all")
        bed_col_dict[species].records.set_index("scaffold", inplace=True)

        bed_col_dict[species].records["color"] = bed_col_dict[species].records["query"].replace(color_df_dict[species]["color"])
        bed_col_dict[species].records.sort_values(by=["scaffold", "start", "end", ], inplace=True)

elif args.synteny_format in ["bed", "bed_with_color"]:
    bed_file_dict = {query: expand_path(bed, skip=not args.expand_paths) for query, bed in zip(query_list, args.input)}

    for species in query_list:
        bed_col_dict[species] = CollectionBED(in_file=get_filenames_for_extension(data_dir_path / species,
                                                                                  extension_list=["bed", "bed.gz"]), header_in_file=True,
                                              format="bed_synteny_track", parsing_mode="all",
                                              scaffold_syn_dict=syn_df_dict[reference] if syn_df_dict[reference] is not None else None,
                                              rename_dict={"query": syn_df_dict[species]} if syn_df_dict[species] is not None else None)

        bed_col_dict[species].records.sort_values(by=["scaffold", "start", "end", ], inplace=True)
        if args.synteny_format != "bed_with_color":
            bed_col_dict[species].records["color"] = bed_col_dict[species].records["query"].replace(color_df_dict[species]["color"])
else:
    raise ValueError("ERROR!!! Unrecognized format of the input file(s)!")

query_scaffold_id_column_name = "query"
query_start_column_name = "query_start"
query_end_column_name = "query_end"
strand_column_name = "strand"
inverted_scaffold_label = "'"


for species in query_list:
    print("Inverting (if necessary) {0} scaffolds...".format(species))
    print("Inverting query coordinates in synteny file...")

    bed_col_dict[species].records = invert_coordinates_in_synteny_table(bed_col_dict[species].records,
                                                                        invertlist_series_dict[species],
                                                                        lenlist_df_dict[species],
                                                                        query_scaffold_id_column_name,
                                                                        query_start_column_name,
                                                                        query_end_column_name,
                                                                        strand_column_name,
                                                                        inverted_scaffold_label)

    lenlist_df_dict[species].rename(index=dict(zip(invertlist_series_dict[species], [scaf + inverted_scaffold_label for scaf in invertlist_series_dict[species]])),
                                    inplace=True)
    #chr_color_df_dict[species].rename(index=dict(zip(invert_list_dict[species], [scaf + label  for scaf in invert_list_dict[species]])),inplace=True)
    orderlist_series_dict[species].replace(dict(zip(invertlist_series_dict[species], [scaf + inverted_scaffold_label for scaf in invertlist_series_dict[species]])),
                                           inplace=True)
    color_df_dict[species].rename(index=dict(zip(invertlist_series_dict[species], [scaf + inverted_scaffold_label for scaf in invertlist_series_dict[species]])),
                                  inplace=True)
#--------------------------------------------------------------------------------------


query_species_color_df_dict = {sp: color_df_dict[sp] for sp in query_list}

print("Lengths of reference chromosomes:")
print(lenlist_df_dict[reference])


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
                                lenlist_df_dict[reference],#reference_scaffold_length_df,
                                orderlist_series_dict[reference],
                                "{0}.initial_min_block_len_{1}".format(args.output_prefix, min_block_length),
                                legend=None if args.hide_legend else Visualization.chromosome_legend(query_species_color_df_dict,
                                                                       orderlist_series_dict[reference]),
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
                                xmax_multiplier=1.3, ymax_multiplier=args.ymax_multiplier,
                                figure_header_height=args.figure_header_height,
                                stranded_tracks=args.stranded,
                                rounded_tracks=args.rounded,
                                stranded_end_tracks=args.stranded_end,
                                xtick_fontsize=args.x_tick_fontsize,
                                subplot_title_fontsize=args.title_fontsize,
                                subplot_title_fontweight='bold'
                                )
    bed_dict_to_xlsx(prefiltered_bed_col_dict, '{0}.initial_min_block_len_{1}'.format(args.output_prefix, min_block_length))

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
                                        lenlist_df_dict[reference],
                                        orderlist_series_dict[reference],
                                        "{0}.{1}".format(args.output_prefix,second_stage_output_suffix),
                                        legend=None if args.hide_legend else Visualization.chromosome_legend(query_species_color_df_dict,
                                                                                                             orderlist_series_dict[reference]),
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
                                        xmax_multiplier=1.3, ymax_multiplier=args.ymax_multiplier,
                                        figure_header_height=args.figure_header_height,
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
                                            lenlist_df_dict[reference],
                                            orderlist_series_dict[reference],
                                            "{0}.{1}.{2}".format(args.output_prefix,
                                                                 second_stage_output_suffix,
                                                                 third_stage_output_suffix),
                                            legend=None if args.hide_legend else Visualization.chromosome_legend(query_species_color_df_dict,
                                                                                                                 orderlist_series_dict[reference]),
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
                                            xmax_multiplier=1.3, ymax_multiplier=args.ymax_multiplier,
                                            stranded_tracks=False,
                                            rounded_tracks=args.rounded,
                                            stranded_end_tracks=False,
                                            xtick_fontsize=args.x_tick_fontsize,
                                            subplot_title_fontsize=args.title_fontsize,
                                            subplot_title_fontweight='bold',
                                            figure_header_height=args.figure_header_height
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
                                                lenlist_df_dict[reference],
                                                orderlist_series_dict[reference],
                                                "{0}.{1}.{2}.{3}".format(args.output_prefix,
                                                                         second_stage_output_suffix,
                                                                         third_stage_output_suffix,
                                                                         forth_stage_output_suffix),
                                                legend=None if args.hide_legend else Visualization.chromosome_legend(query_species_color_df_dict,
                                                                                                                     orderlist_series_dict[reference]),
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
                                                xmax_multiplier=1.3, ymax_multiplier=args.ymax_multiplier,
                                                figure_header_height=args.figure_header_height,
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
