#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os

import argparse
import glob
from functools import partial
from pathlib import Path

import pandas as pd
import numpy as np
from copy import deepcopy
from distinctipy import distinctipy

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

from MACE.Visualization.Polygons import LinearChromosome
from MACE.Visualization.Connectors import CubicBezierConnector
from RouToolPa.Parsers.PSL import CollectionPSL


def split_comma_separated_list(string):
    return string.split(",")


def rgb_tuple_to_hex(rgb_tuple):
    color_code = "#"
    for i in [0, 1, 2]:
        color_code += "{:02X}".format(int(255 * rgb_tuple[i]))

    return color_code


def get_filenames_for_extension(dir_path, extension_list, force_uniq=True):
    filelist = []
    for extension in extension_list:
        filelist += list(glob.glob(str(dir_path) + "/*{0}".format(extension)))
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
    columns_list = list(temp_df.columns)
    temp_df["length_column"] = 0
    temp_df.set_index(scaffold_column, inplace=True)
    #print (temp_df)
    #temp_df.to_csv("tmp", sep="\t", index=True, header=True)
    #print(scaffold_list)
    for scaffold in temp_df.index.unique():
        #print(scaffold)
        #print(temp_df)
        #print(length_df)
        temp_df.loc[scaffold, "length_column"] = length_df.loc[scaffold, "length"]

    temp_df.loc[temp_df.index.isin(scaffold_list), start_column], temp_df.loc[temp_df.index.isin(scaffold_list), end_column] = temp_df.loc[temp_df.index.isin(scaffold_list), "length_column"] - temp_df.loc[temp_df.index.isin(scaffold_list), end_column], \
               temp_df.loc[temp_df.index.isin(scaffold_list), "length_column"] - temp_df.loc[temp_df.index.isin(scaffold_list), start_column]

    plus_indexes, minus_indexes = (temp_df[strand_column] == "+") & temp_df.index.isin(scaffold_list), (temp_df[strand_column] == "-") & temp_df.index.isin(scaffold_list)
    temp_df.loc[plus_indexes, strand_column], temp_df.loc[minus_indexes, strand_column] = "-", "+"
    temp_df.reset_index(drop=False, inplace=True)
    if inverted_scaffolds_label is not None:
        for scaffold in scaffold_list:
            temp_df.loc[temp_df[scaffold_column] == scaffold, scaffold_column] = scaffold + inverted_scaffolds_label
    return temp_df[columns_list]  # remove added length column and restore column order


def invert_coordinates_in_region_table(df, scaffold_list, length_df, scaffold_column, start_column, end_column, inverted_scaffolds_label="'"):
    if df.empty:
        return df
    temp_df = deepcopy(df)
    #print(df)
    #print(length_df)
    if temp_df.index.name != scaffold_column:
        temp_df.reset_index(inplace=True, drop=False)
        temp_df.set_index(scaffold_column, inplace=True)

    columns_list = list(temp_df.columns)

    for scaffold in temp_df.index.unique():
        #print(scaffold)
        #print(temp_df)
        #print(length_df)
        #print(length_df)
        temp_df.loc[scaffold, "length_column"] = length_df.loc[scaffold, "length"]

    temp_df.loc[temp_df.index.isin(scaffold_list), start_column], temp_df.loc[temp_df.index.isin(scaffold_list), end_column] = temp_df.loc[temp_df.index.isin(scaffold_list), "length_column"] - temp_df.loc[temp_df.index.isin(scaffold_list), end_column], \
               temp_df.loc[temp_df.index.isin(scaffold_list), "length_column"] - temp_df.loc[temp_df.index.isin(scaffold_list), start_column]

    temp_df.reset_index(drop=False, inplace=True)
    if inverted_scaffolds_label is not None:
        for scaffold in scaffold_list:
            temp_df.loc[temp_df[scaffold_column] == scaffold, scaffold_column] = scaffold + inverted_scaffolds_label

    #print(temp_df)
    temp_df.set_index(scaffold_column, inplace=True)
    return temp_df[columns_list]  # remove added length column and restore column order


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True,
                    help="Input directory with data. Must contain subfolder for each genome. "
                         "Subfolders should have same name as genomes in --genome_orderlist"
                         "Each subfolder should contain: *.whitelist, *.len and synteny file "
                         "(except for the last genome). *.orderlist, *.invertlist and *.syn file are optional")
parser.add_argument("--genome_orderlist", action="store", dest="genome_orderlist", required=True,
                    type=split_comma_separated_list,
                    help="Comma-separated list of genomes to be used in figure.")
parser.add_argument("--genome_labellist", action="store", dest="genome_labellist", default=None,
                    type=split_comma_separated_list,
                    help="Comma-separated list of genome labels to be used in figure instead of genome names. "
                         "Must follow the same order as --genome_orderlist. If not set genome names will serve as labels.")
parser.add_argument("--invert_genome_order", action="store_true", dest="invert_genome_order", default=False,
                    help="Invert order of the genomes in the --genome_orderlist. Default: False")
parser.add_argument("--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in synonym file. Default: 0")
parser.add_argument("--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in synonym file synonym. Default: 1")
parser.add_argument("--synteny_format", action="store", dest="synteny_format", default="psl",
                    help="Format of the synteny file. Allowed: psl(default)")
parser.add_argument("--remove_nested_blocks", action="store_true", dest="remove_nested_blocks", default=False,
                    help="Remove nested blocks. To be removed block must be nested in other block either on target or query side. Default: False")
parser.add_argument("--remove_same_coords_blocks", action="store_true", dest="remove_same_coords_blocks", default=False,
                    help="Remove blocks with exactly the same coordinates on query or target side. Default: False")
parser.add_argument("--min_len_threshold", action="store", dest="min_len_threshold", default=0, type=int,
                    help="Minimum length of rearranged block to be highlighted. "
                         "Recommended value for mammalian-size genomes ranges between 200'000 and 1000'000"
                         "Default: 0, i.e. all rearranged blocks will be highlighted")
parser.add_argument("--scaffold_prefix_cut", action="store", dest="scaffold_prefix_cut", default=3, type=int,
                    help="Length of prefix to be cut from every scaffold id. "
                         "Default: 3, i.e 'chr3' will be cut to '3' on the figure. Set this option to zero to avoid cutting.")
parser.add_argument("--invert_major_strand", action="store_true", dest="invert_major_strand", default=False,
                    help="Invert major strand, i.e treat all inverted blocks as normal, and vice versa. "
                         "This flag affects all genomes and all chromosomes."
                         "If you wish to invert the major strand for the specific genome or even for a specific "
                         "chromosome, please, use genome=specific swithstrandlist file.  Default: False")

parser.add_argument("--inversion_color", action="store", dest="inversion_color", default="red",
                    help="Color to use to highlight inversions on the plot. Must be a color recognized by Matplotlib. Default: 'red'")
parser.add_argument("--translocation_color", action="store", dest="translocation_color", default="blue",
                    help="Color to use to highlight translocation on the plot. Must be a color recognized by Matplotlib. Default: 'blue'")
parser.add_argument("--default_color", action="store", dest="default_color", default='default',
                    help="Color to use for connectors between synteny blocks on the plot. Must be a color recognized by Matplotlib. Default: 'lightgrey'")

parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=split_comma_separated_list,
                    default=("png", "svg"),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg,png")

parser.add_argument("-l", "--title", action="store", dest="title", default="Macrosynteny",
                    help="Suptitle of figure. Default: 'Macrosynteny'")
parser.add_argument("--title_fontsize", action="store", dest="title_fontsize", default=20, type=int,
                    help="Fontsize of the figure. Default: 20")
parser.add_argument("--hide_chromosome_labels", action="store_true", dest="hide_chromosome_labels", default=False,
                    help="Hide chromosome labels. Default: False")
parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float, default=0.05,
                    help="Adjust left border of subplots on the figure. Default: matplotlib defaults")

parser.add_argument("--subplots_adjust_right", action="store", dest="subplots_adjust_right", type=float, default=0.98,
                    help="Adjust right border of subplots on the figure. Default: matplotlib defaults")

parser.add_argument("--subplots_adjust_top", action="store", dest="subplots_adjust_top", type=float, default=0.90,
                    help="Adjust top border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_bottom", action="store", dest="subplots_adjust_bottom", type=float, default=0.05,
                    help="Adjust bottom border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--figure_width", action="store", dest="figure_width", type=float, default=8,
                    help="Width of figure in inches. Default: 8")
parser.add_argument("--figure_height_per_genome", action="store", dest="figure_height_per_genome",
                    type=float, default=1,
                    help="Height of figure per genome track. Default: 1")
parser.add_argument("--chromosome_label_fontsize", action="store", dest="chromosome_label_fontsize", type=int, default=7,
                    help="Fontsize of chromosome labels. Default: 7")
parser.add_argument("--genome_label_fontsize", action="store", dest="genome_label_fontsize", type=int, default=11,
                    help="Fontsize of genome labels. Default: 11")
parser.add_argument("--chromosome_label_angle", action="store", dest="chromosome_label_angle", type=int, default=0,
                    help="Angle to rotate chromosome labels on the plot. Default: 0")
parser.add_argument("--genome_label_angle", action="store", dest="genome_label_angle", type=int, default=0,
                    help="Angle to rotate genome labels on the plot. Default: 0")
parser.add_argument("--genome_distance", action="store", dest="genome_distance", type=int, default=100,
                    help="Distance between genomes on the plot. Default: 100")
parser.add_argument("--smooth_multiplicator", action="store", dest="smooth_multiplicator", type=float, default=4,
                    help="Multiplicator used to control smoothing of the chromosome ends and visual width of centromere. "
                         "Reduction of it will narrow centromere and make chromosomes closer to the rectangle."
                         "Default value (4) is good for usual carnivora genomes (2-3 Gbp, with 2n~16-50)")


#parser.add_argument("--subplot_scale", action="store_true", dest="subplot_scale",
#                    help="Scale feature x size by subplot x/y ratio. Default: off")
#parser.add_argument("--track_group_scale", action="store_true", dest="track_group_scale",
#                    help="Scale feature x size by track_group x/y ratio. Default: off")
#
#parser.add_argument("--x_tick_fontsize", action="store", dest="x_tick_fontsize", type=int, default=None,
#                    help="Fontsize of xticks. Default: matplotlib default")
#parser.add_argument("--invert_coordinates_for_target_negative_strand", action="store_true",
#                    dest="invert_coordinates_for_target_negative_strand",
#                    default=False,
#                    help="Invert coordinates for target negative strand. Default: False")


args = parser.parse_args()

data_dir = args.input_dir
data_dir_path = Path(data_dir)

genome_orderlist = args.genome_orderlist[::-1] if args.invert_genome_order else args.genome_orderlist
if args.genome_labellist is None:
    genome_labellist = genome_orderlist
else:
    genome_labellist = args.genome_labellist[::-1] if args.invert_genome_order else args.genome_labellist

syn_file_key_column, syn_file_value_column = args.syn_file_key_column, args.syn_file_value_column

#label_dict = {genome: label for genome, label in zip(args.genome_orderlist,
#                                                     args.genome_orderlist if args.genome_labellist is None else args.genome_labellist)}

synteny_format = args.synteny_format

inverted_scaffold_label = "'"

# read files
print(data_dir_path)
# whitelist is obligatory
for genome in genome_orderlist:
    print(data_dir_path / genome)
    #print(get_filenames_for_extension(data_dir_path / genome, extension_list=["whitelist"]))

whitelist_series_dict = {}
orderlist_series_dict = {}
invertlist_series_dict = {}
lenlist_df_dict = {}
syn_df_dict = {}
color_df_dict = {}
queryswithstrandlist_series_dict = {}
targetswithstrandlist_series_dict = {}
centromere_df_dict = {}
for genome in genome_orderlist:
    # nonempty whitelist file is necessary for each genome
    try:
        whitelist_series_dict[genome] = pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["whitelist"]),
                                                    sep="\t", header=None, comment="#").squeeze("columns")
    except pd.errors.EmptyDataError:
        raise pd.errors.EmptyDataError("ERROR!!! Whitelist for {0} is empty. Add relevant scaffold ids to it!")

    # nonempty lenlist file is necessary for each genome
    try:
        lenlist_df_dict[genome] = pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["len"]),
                                             sep="\t", header=None, comment="#", names=["scaffold", "length"], index_col=0)
    except pd.errors.EmptyDataError:
        raise pd.errors.EmptyDataError("ERROR!!! lenlist for {0} is empty. add scaffold ids and its lengths to it!")

    # orderlist might be empty or absent
    if get_filenames_for_extension(data_dir_path / genome, extension_list=["orderlist"]) is None:
        orderlist_series_dict[genome] = pd.Series(dtype=str)
    else:
        try:
            orderlist_series_dict[genome] = pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["orderlist"]),
                                                        sep="\t", header=None, comment="#").squeeze("columns")
        except pd.errors.EmptyDataError:
            orderlist_series_dict[genome] = pd.Series(dtype=str)

    # invertlist might be empty or absent
    if get_filenames_for_extension(data_dir_path / genome, extension_list=["invertlist"]) is None:
        invertlist_series_dict[genome] = pd.Series(dtype=str)
    else:
        try:
            invertlist_series_dict[genome] = pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["invertlist"]),
                                                         sep="\t", header=None, comment="#").squeeze("columns")
        except pd.errors.EmptyDataError:
            invertlist_series_dict[genome] = pd.Series(dtype=str)

    # queryswitchstrandlist might be empty or absent
    if get_filenames_for_extension(data_dir_path / genome, extension_list=["queryswitchstrandlist"]) is None:
        queryswithstrandlist_series_dict[genome] = pd.Series(dtype=str)
    else:
        try:
            queryswithstrandlist_series_dict[genome] = pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["queryswitchstrandlist"]),
                                                                   sep="\t", header=None, comment="#").squeeze("columns")
        except pd.errors.EmptyDataError:
            queryswithstrandlist_series_dict[genome] = pd.Series(dtype=str)

    # targetswitchstrandlist might be empty or absent
    if get_filenames_for_extension(data_dir_path / genome, extension_list=["targetswitchstrandlist"]) is None:
        targetswithstrandlist_series_dict[genome] = pd.Series(dtype=str)
    else:
        try:
            targetswithstrandlist_series_dict[genome] = pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["targetswitchstrandlist"]),
                                                                    sep="\t", header=None, comment="#").squeeze("columns")
        except pd.errors.EmptyDataError:
            targetswithstrandlist_series_dict[genome] = pd.Series(dtype=str)

    # syn file might be empty or absent
    if get_filenames_for_extension(data_dir_path / genome, extension_list=["syn"]) is None:
        syn_df_dict[genome] = pd.DataFrame(columns=["syn"])
    else:
        try:
            syn_df_dict[genome] = pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["syn"]), usecols=(syn_file_key_column, syn_file_value_column),
                                              sep="\t", header=None, comment="#",
                                              names=["key", "syn"] if syn_file_key_column <= syn_file_value_column else ["syn", "key"]).set_index("key")
        except pd.errors.EmptyDataError:
            syn_df_dict[genome] = pd.DataFrame(columns=["syn"])

    # centromere.bed might be empty or absent
    if get_filenames_for_extension(data_dir_path / genome, extension_list=["centromere.bed"]) is None:
        centromere_df_dict[genome] = pd.DataFrame(columns=["scaffold_id", "start", "end"])
        centromere_df_dict[genome].set_index("scaffold_id", inplace=True)
    else:
        try:
            centromere_df_dict[genome] = pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["centromere.bed"]),
                                                     usecols=(0, 1, 2),
                                                     index_col=0,
                                                     header=None,
                                                     sep="\t", names=["scaffold_id", "start", "end"])
        except pd.errors.EmptyDataError:
            centromere_df_dict[genome] = pd.DataFrame(columns=["scaffold_id", "start", "end"])
            centromere_df_dict[genome].set_index("scaffold_id", inplace=True)
    centromere_df_dict[genome].rename(index=syn_df_dict[genome]["syn"].to_dict(), inplace=True)

    # color file might be empty or absent
    if get_filenames_for_extension(data_dir_path / genome, extension_list=["colors"]) is None:
        color_df_dict[genome] = pd.DataFrame(columns=["color"])
    else:
        try:
            color_df_dict[genome] = pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["colors"]),
                                             sep="\t", header=0, names=["scaffold", "color"], index_col=0)
        except pd.errors.EmptyDataError:
            color_df_dict[genome] = pd.DataFrame(columns=["scaffold_id", "color"])
            color_df_dict[genome].set_index("scaffold_id", inplace=True)


#whitelist_series_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["whitelist"]),
#                                             sep="\t", header=None, comment="#").squeeze("columns") for genome in genome_orderlist}
# orderlist might be absent in folders
#orderlist_series_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["orderlist"]),
#                                             sep="\t", header=None, comment="#").squeeze("columns") if get_filenames_for_extension(data_dir_path / genome, extension_list=["orderlist"]) is not None else pd.Series(dtype=str)  for genome in genome_orderlist}
# invertlist might be absent in folders
#invertlist_series_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["invertlist"]),
#                                             sep="\t", header=None, comment="#").squeeze("columns") if get_filenames_for_extension(data_dir_path / genome, extension_list=["invertlist"]) is not None else pd.Series(dtype=str) for genome in genome_orderlist}
# lenlist is obligatory
#lenlist_df_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["len"]),
#                                             sep="\t", header=None, comment="#", names=["scaffold", "length"], index_col=0) for genome in genome_orderlist}
# synfile might be absent
#syn_df_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["syn"]), usecols=(syn_file_key_column, syn_file_value_column),
#                                             sep="\t", header=None, comment="#",
#                                   names=["key", "syn"] if syn_file_key_column <= syn_file_value_column else ["syn", "key"]).set_index("key") if get_filenames_for_extension(data_dir_path / genome, extension_list=["syn"]) is not None else pd.DataFrame(columns=["syn"]) for genome in genome_orderlist}

# swithstrandlist might be absent in folders
#swithstrandlist_series_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["switchstrandlist"]),
#                                             sep="\t", header=None, comment="#").squeeze("columns") if get_filenames_for_extension(data_dir_path / genome, extension_list=["switchstrandlist"]) is not None else pd.Series(dtype=str) for genome in genome_orderlist}


#print("AAAAAAAAAA")
#print(lenlist_df_dict)
#filter len list
for genome in genome_orderlist:
    lenlist_df_dict[genome] = lenlist_df_dict[genome].loc[lenlist_df_dict[genome].index.isin(whitelist_series_dict[genome])]
#print("BBBBBBBBBB")
#print(lenlist_df_dict)
# rename len list according to synonyms
for genome in genome_orderlist:
    if not syn_df_dict[genome].empty:
        lenlist_df_dict[genome].rename(index=syn_df_dict[genome]["syn"].to_dict(), inplace=True)
#print("CCCCCCCCCCCCCCCCCCCCC")
# reorder len list according the orderlist
#print(lenlist_df_dict)
for genome in genome_orderlist:
    lenlist_df_dict[genome] = lenlist_df_dict[genome].reindex(orderlist_series_dict[genome]).dropna()
#print("DDDDDDDDDDDDDDDD")
#print(lenlist_df_dict)
# count number of chromosomes/scaffold and total length of them
total_len_dict = {genome: sum(lenlist_df_dict[genome]["length"]) for genome in genome_orderlist}
chr_number_dict = {genome: len(lenlist_df_dict[genome]) for genome in genome_orderlist}

max_genome_length = max(list(total_len_dict.values()))
#print("CCCCCCCCCCC")
#print(queryswithstrandlist_series_dict)

if synteny_format == "psl":
    # psl or psl.gz files must exist for all genomes except the last one
    strand_column_name = "strand"

    query_scaffold_id_column_name = "qName"
    query_start_column_name = "qStart"
    query_end_column_name = "qEnd"
    query_block_len_column_name = "qHitLen"

    target_scaffold_id_column_name = "tName"
    target_start_column_name = "tStart"
    target_end_column_name = "tEnd"
    target_block_len_column_name = "tHitLen"

    connector_color_idx = None
    strand_idx = 0

    target_scaffold_idx = 4
    target_start_idx = 5
    target_end_idx = 6

    query_scaffold_idx = 1
    query_start_idx = 2
    query_end_idx = 3

    synteny_dict = {}
    for genome_index in range(0, len(genome_orderlist)-1):
        print("\n")
        print("Filename: {0}".format(get_filenames_for_extension(data_dir_path / genome_orderlist[genome_index],
                                                                                                  extension_list=["psl", "psl.gz"])))
        print("Query: {0}".format(genome_orderlist[genome_index]))
        print("Target: {0}".format(genome_orderlist[genome_index + 1]))
        synteny_dict[genome_orderlist[genome_index]] = CollectionPSL(get_filenames_for_extension(data_dir_path / genome_orderlist[genome_index],
                                                                                                  extension_list=["psl", "psl.gz"]),
                                                                     target_white_list=whitelist_series_dict[genome_orderlist[genome_index + 1]],
                                                                     query_white_list=whitelist_series_dict[genome_orderlist[genome_index]],
                                                                     #target_syn_dict=syn_df_dict[genome_orderlist[genome_index + 1]]["syn"].to_dict(),
                                                                     #query_syn_dict=syn_df_dict[genome_orderlist[genome_index]]["syn"].to_dict(),
                                                                     parsing_mode="coordinates_only").records.sort_values(by=[query_scaffold_id_column_name,
                                                                                                                              query_start_column_name,
                                                                                                                              query_end_column_name,
                                                                                                                              target_scaffold_id_column_name,
                                                                                                                              target_start_column_name,
                                                                                                                              target_end_column_name])
else:
    raise ValueError("ERROR!!! {0} format is not implemented yet!".format("psl"))

#----------------------------- Assignment of IDs to synteny blocks -------------------------------
for genome in synteny_dict:
    synteny_dict[genome]['synteny_block_id'] = ["SB_{0}".format(block_id) for block_id in range(1, len(synteny_dict[genome]) + 1)]
    #print(synteny_dict[genome])
#-------------------------------------------------------------------------------------------------

#print(syn_df_dict)
#print(synteny_dict["dingo"])
print("\n\n")
genome_number = len(genome_orderlist)

#------------------ Preprocessing of synteny blocks -----------------------------------
##----------------------- Detect nested blocks ----------------------------------------


def detect_nested_blocks(df,
                         nested_in_block_column_name,
                         #query_same_cooords_in_block_column_name,
                         query_nested_in_block_column_name,
                         query_scaffold_id_column_name,
                         query_start_column_name, query_end_column_name,
                         #target_same_cooords_in_block_column_name,
                         target_nested_in_block_column_name,
                         target_scaffold_id_column_name,
                         target_start_column_name, target_end_column_name):
    sorted_df = df.sort_values(by=[query_scaffold_id_column_name, query_start_column_name, query_end_column_name])
    for column_name in query_nested_in_block_column_name, target_nested_in_block_column_name: #, query_same_cooords_in_block_column_name, target_same_cooords_in_block_column_name:
        sorted_df[column_name] = pd.NA
    #sorted_df[query_nested_in_block_column_name] = pd.NA
    #sorted_df[target_nested_in_block_column_name] = pd.NA

    for row_index in range(0, len(sorted_df)):
        block_start = sorted_df.iloc[row_index][query_start_column_name]
        block_end = sorted_df.iloc[row_index][query_end_column_name]
        #check_df = sorted_df.iloc[:row_index]
        #print(sorted_df.iloc[row_index])
        #print(check_df)
        nested_in_block_set = set(sorted_df.iloc[:row_index][sorted_df.iloc[:row_index][query_end_column_name] >= block_end]['synteny_block_id'])
        nested_in_block_set |= set(sorted_df.iloc[row_index+1:][(sorted_df.iloc[row_index+1:][query_start_column_name] == block_start) & (sorted_df.iloc[row_index+1:][query_end_column_name] >= block_end)]['synteny_block_id'])

        #print(nested_in_block_set)
        if nested_in_block_set:
            sorted_df[query_nested_in_block_column_name].iloc[row_index] = ",".join(nested_in_block_set)

        #same_coords_set = set(sorted_df[(sorted_df[query_start_column_name] == block_start) & (sorted_df[query_end_column_name] == block_end)]['synteny_block_id'])
        #same_coords_set.remove(sorted_df.iloc[row_index]['synteny_block_id'])

        #if same_coords_set:
        #    sorted_df[query_same_cooords_in_block_column_name].iloc[row_index] = ",".join(same_coords_set)

    sorted_df = sorted_df.sort_values(by=[target_scaffold_id_column_name, target_start_column_name, target_end_column_name])
    for row_index in range(0, len(sorted_df)):
        block_start = sorted_df.iloc[row_index][target_start_column_name]
        block_end = sorted_df.iloc[row_index][target_end_column_name]
        #check_df = sorted_df.iloc[:row_index]
        nested_in_block_set = set(sorted_df.iloc[:row_index][sorted_df.iloc[:row_index][target_end_column_name] >= block_end]['synteny_block_id'])
        nested_in_block_set |= set(sorted_df.iloc[row_index+1:][(sorted_df.iloc[row_index+1:][target_start_column_name] == block_start) & (sorted_df.iloc[row_index+1:][target_end_column_name] >= block_end)]['synteny_block_id'])

        if nested_in_block_set:
            sorted_df[target_nested_in_block_column_name].iloc[row_index] = ",".join(nested_in_block_set)

        #same_coords_set = set(sorted_df[(sorted_df[target_start_column_name] == block_start) & (sorted_df[target_end_column_name] == block_end)]['synteny_block_id'])
        #same_coords_set.remove(sorted_df.iloc[row_index]['synteny_block_id'])
        #if same_coords_set:
        #    sorted_df[target_same_cooords_in_block_column_name].iloc[row_index] = ",".join(same_coords_set)

    def get_nested(row):
        if row.hasnans:
            return pd.NA
        else:
            return ",".join(set(row.iloc[0].split(",")) & set(row.iloc[1].split(",")))
    sorted_df[nested_in_block_column_name] = sorted_df[[query_nested_in_block_column_name, target_nested_in_block_column_name]].apply(get_nested, axis=1)
    #print(sorted_df)
    return sorted_df


def detect_same_coords_blocks(df,
                              query_same_coords_in_block_column_name,
                              query_scaffold_id_column_name,
                              query_start_column_name, query_end_column_name,
                              target_same_coords_in_block_column_name,
                              target_scaffold_id_column_name,
                              target_start_column_name, target_end_column_name):
    output_df = deepcopy(df)
    for column_name in query_same_coords_in_block_column_name, target_same_coords_in_block_column_name:
        output_df[column_name] = pd.NA
    #print(output_df)
    for row_index in range(0, len(output_df)):
        query_same_coords_set = set(output_df[(output_df[query_scaffold_id_column_name] == output_df.iloc[row_index][query_scaffold_id_column_name]) & \
                                              (output_df[query_start_column_name] == output_df.iloc[row_index][query_start_column_name]) & \
                                              (output_df[query_end_column_name] == output_df.iloc[row_index][query_end_column_name])]['synteny_block_id'])
        query_same_coords_set.remove(output_df.iloc[row_index]['synteny_block_id'])

        if query_same_coords_set:
            output_df[query_same_coords_in_block_column_name].iloc[row_index] = ",".join(query_same_coords_set)
            
        target_same_coords_set = set(output_df[(output_df[target_scaffold_id_column_name] == output_df.iloc[row_index][target_scaffold_id_column_name]) & \
                                              (output_df[target_start_column_name] == output_df.iloc[row_index][target_start_column_name]) & \
                                              (output_df[target_end_column_name] == output_df.iloc[row_index][target_end_column_name])]['synteny_block_id'])
        target_same_coords_set.remove(output_df.iloc[row_index]['synteny_block_id'])

        if target_same_coords_set:
            output_df[target_same_coords_in_block_column_name].iloc[row_index] = ",".join(target_same_coords_set)

    return output_df


def detect_overlapping_blocks(df,
                              query_overlapping_block_column_name,
                              query_overlapping_fraction_column_name,
                              query_reverse_overlapping_fraction_column_name,
                              query_start_column_name, query_end_column_name,
                              target_overlapping_block_column_name,
                              target_overlapping_fraction_column_name,
                              target_reverse_overlapping_fraction_column_name,
                              target_start_column_name, target_end_column_name):
    output_df = deepcopy(df)
    for column_name in query_overlapping_block_column_name, target_overlapping_block_column_name:
        output_df[column_name] = pd.NA
    for column_name in query_overlapping_fraction_column_name, \
                       query_reverse_overlapping_fraction_column_name, \
                       target_overlapping_fraction_column_name, \
                       target_reverse_overlapping_fraction_column_name:
        output_df[column_name] = np.nan

    #print(output_df)
    for row_index in range(0, len(output_df)):
        output_df["query.t_e - x_s"] = output_df[query_end_column_name].iloc[row_index] - output_df[query_start_column_name]

        output_df["query.t_e - x_s"] = output_df[query_end_column_name].iloc[row_index] - output_df[query_start_column_name]

        output_df["query.t_e - x_s"] = output_df[query_end_column_name].iloc[row_index] - output_df[query_start_column_name]

        output_df["query.t_e - x_s"] = output_df[query_end_column_name].iloc[row_index] - output_df[query_start_column_name]

    return output_df


detect_nested_blocks_preset = partial(detect_nested_blocks,
                                      nested_in_block_column_name="nested_in",
                                      query_nested_in_block_column_name="query_nested_in",
                                      #query_same_cooords_in_block_column_name="query_same_coords",
                                      query_scaffold_id_column_name=query_scaffold_id_column_name,
                                      query_start_column_name=query_start_column_name,
                                      query_end_column_name=query_end_column_name,
                                      target_nested_in_block_column_name="target_nested_in",
                                      #target_same_cooords_in_block_column_name="target_same_coords",
                                      target_scaffold_id_column_name=target_scaffold_id_column_name,
                                      target_start_column_name=target_start_column_name,
                                      target_end_column_name=target_end_column_name)

tmp_dict = {}
block_remove_dict = {}
for genome_index in range(0, genome_number - 1):
    genome = genome_orderlist[genome_index]
    block_remove_dict[genome] = []
    tmp_dict[genome] = detect_same_coords_blocks(synteny_dict[genome],
                                                 query_same_coords_in_block_column_name="query_same_coords",
                                                 query_scaffold_id_column_name=query_scaffold_id_column_name,
                                                 query_start_column_name=query_start_column_name,
                                                 query_end_column_name=query_end_column_name,
                                                 target_same_coords_in_block_column_name="target_same_coords",
                                                 target_scaffold_id_column_name=target_scaffold_id_column_name,
                                                 target_start_column_name=target_start_column_name,
                                                 target_end_column_name=target_end_column_name)
    tmp_dict[genome] = tmp_dict[genome].groupby(by=[query_scaffold_id_column_name, target_scaffold_id_column_name],
                                                sort=False, group_keys=False).apply(detect_nested_blocks_preset)

    #tmp_dict[genome].sort_values(by=[query_scaffold_id_column_name,
    #                                 query_start_column_name,
    #                                 query_end_column_name,
    #                                 target_scaffold_id_column_name,
    #                                 target_start_column_name,
    #                                 target_end_column_name]).to_csv("TMP_{0}.{1}.to.{2}.raw.tab".format(args.output_prefix,
    #                                                                 genome,
    #                                                                 genome_orderlist[genome_index + 1]),
    #                              sep="\t", index=False, header=True)
    if args.remove_same_coords_blocks:
        block_remove_dict[genome] = list(tmp_dict[genome]["target_same_coords"].dropna())

    if args.remove_nested_blocks:
        block_remove_dict[genome] += list(tmp_dict[genome][tmp_dict[genome]["query_nested_in"].notna()]["synteny_block_id"])
        block_remove_dict[genome] += list(tmp_dict[genome][tmp_dict[genome]["target_nested_in"].notna()]["synteny_block_id"])

    block_remove_dict[genome] = set(block_remove_dict[genome])
    synteny_dict[genome] = synteny_dict[genome][~synteny_dict[genome]["synteny_block_id"].isin(block_remove_dict[genome])]


#print(block_remove_dict)
#------------------ Classification of synteny blocks ----------------------------------
for genome_index in range(0, genome_number - 1): # genome_orderlist[:-1]:  # all query genomes
    target_genome_index = genome_index + 1
    genome = genome_orderlist[genome_index]
    target_genome = genome_orderlist[target_genome_index]
    #pd.set_option('display.max_rows', 40)
    columns_list = list(synteny_dict[genome].columns)
    hit_sum = synteny_dict[genome][[strand_column_name,
                                    query_scaffold_id_column_name,
                                    target_scaffold_id_column_name,
                                    query_block_len_column_name]].groupby(by=[query_scaffold_id_column_name,
                                                                              target_scaffold_id_column_name,
                                                                              strand_column_name]).sum()

    #print(synteny_dict[genome])
    synteny_dict[genome]["type"] = "normal"
    synteny_dict[genome]["connector_color"] = args.default_color
    synteny_dict[genome]["connector_zorder"] = 0
    synteny_dict[genome].set_index([query_scaffold_id_column_name, target_scaffold_id_column_name], inplace=True)
    #print(synteny_dict[genome])
    # major_strand_series is a series with two level index(qName, tName)
    major_strand_series = hit_sum.groupby(by=[query_scaffold_id_column_name,
                                              target_scaffold_id_column_name]).apply(lambda df: df.idxmax()[query_block_len_column_name][2])
    #print("AAAA")
    #print(genome)
    #print(major_strand_series)
    #print("BBBB")
    minus_strand_index = major_strand_series == "-"
    plus_strand_index = major_strand_series == "+"

    if args.invert_major_strand:  # first apply a global major strand switch

        print("Flag --invert_major_strand flag was set. Switching major strand for all genomes and all chromosomes...")
        major_strand_series[minus_strand_index] = "+"
        major_strand_series[plus_strand_index] = "-"

    # apply a local major strand switch

    if not queryswithstrandlist_series_dict[genome].empty:
        print(f"Switching major strand for {genome} as query for {genome} vs {target_genome} alignment...")
        #print(major_strand_series)
        #switch_strand_scaffolds = set(swithstrandlist_series_dict[genome]) & set(major_strand_series.index.get_level_values(0))
        #print(major_strand_series == "-")
        switch_series = pd.Series(major_strand_series.index.get_level_values(query_scaffold_id_column_name).isin(queryswithstrandlist_series_dict[genome]))
        switch_series.index = major_strand_series.index

        minus_strand_index = switch_series & (major_strand_series == "-")
        plus_strand_index = switch_series & (major_strand_series == "+")
        major_strand_series[minus_strand_index] = "+"
        major_strand_series[plus_strand_index] = "-"

    # apply switchstrand list of target genome.
    print((genome, target_genome))
    #print(targetswithstrandlist_series_dict[target_genome])
    if not targetswithstrandlist_series_dict[target_genome].empty:
        print(f"Switching major strand for {target_genome} as target for {genome} vs {target_genome} alignment...")
        switch_series = pd.Series(major_strand_series.index.get_level_values(target_scaffold_id_column_name).isin(targetswithstrandlist_series_dict[target_genome]))
        switch_series.index = major_strand_series.index

        minus_strand_index = switch_series & (major_strand_series == "-")
        plus_strand_index = switch_series & (major_strand_series == "+")
        major_strand_series[minus_strand_index] = "+"
        major_strand_series[plus_strand_index] = "-"

    synteny_dict[genome]["major_strand"] = major_strand_series
    #print("CCCC")
    #print(major_strand_series)
    synteny_dict[genome].reset_index(level=1, drop=False, inplace=True)
    #print(genome)
    #print(hit_sum)
    #detect translocations from query side
    major_query_homolog_series = hit_sum.droplevel(level=2).groupby(by=[query_scaffold_id_column_name,
                                                                        target_scaffold_id_column_name]).sum().groupby(by=[query_scaffold_id_column_name]).apply(lambda df: df.idxmax()[query_block_len_column_name][1])

    synteny_dict[genome]["major_query_homolog"] = major_query_homolog_series
    synteny_dict[genome].loc[synteny_dict[genome][target_scaffold_id_column_name] != synteny_dict[genome]["major_query_homolog"], "connector_color"] = args.translocation_color
    synteny_dict[genome].loc[synteny_dict[genome][target_scaffold_id_column_name] != synteny_dict[genome]["major_query_homolog"], "type"] = "translocation"
    synteny_dict[genome].loc[synteny_dict[genome][target_scaffold_id_column_name] != synteny_dict[genome]["major_query_homolog"], "connector_zorder"] = 50

    #detect translocations from target side
    synteny_dict[genome].reset_index(level=0, drop=False, inplace=True)
    #print(hit_sum)
    major_target_homolog_series = hit_sum.droplevel(level=2).groupby(by=[target_scaffold_id_column_name,
                                                                         query_scaffold_id_column_name]).sum().groupby(by=[target_scaffold_id_column_name]).apply(lambda df: df.idxmax()[query_block_len_column_name][1])
    #print(major_target_homolog_series)

    synteny_dict[genome].set_index(target_scaffold_id_column_name, inplace=True)
    synteny_dict[genome]["major_target_homolog"] = major_target_homolog_series
    synteny_dict[genome].loc[synteny_dict[genome][query_scaffold_id_column_name] != synteny_dict[genome]["major_target_homolog"], "connector_color"] = args.translocation_color
    synteny_dict[genome].loc[synteny_dict[genome][query_scaffold_id_column_name] != synteny_dict[genome]["major_target_homolog"], "type"] = "translocation"
    synteny_dict[genome].loc[synteny_dict[genome][query_scaffold_id_column_name] != synteny_dict[genome]["major_target_homolog"], "connector_zorder"] = 50
    #synteny_dict[genome].to_csv("AAAAAAAAAAAAAA.{0}.tmp".format(genome), sep="\t")
    #print(synteny_dict[genome])

    synteny_dict[genome].reset_index(level=0, drop=False, inplace=True)

    synteny_dict[genome].loc[synteny_dict[genome][strand_column_name] != synteny_dict[genome]["major_strand"], "connector_color"] = args.inversion_color
    synteny_dict[genome].loc[synteny_dict[genome][strand_column_name] != synteny_dict[genome]["major_strand"], "type"] = "inversion"
    synteny_dict[genome].loc[synteny_dict[genome][strand_column_name] != synteny_dict[genome]["major_strand"], "connector_zorder"] = 20

    connector_color_idx = list(synteny_dict[genome].columns).index("connector_color")  # len(columns_list) + 1
    connector_zorder_idx = list(synteny_dict[genome].columns).index("connector_zorder")  # len(columns_list) + 2
    synteny_dict[genome] = synteny_dict[genome][columns_list + ["type", "connector_color", "connector_zorder"]]

    if args.min_len_threshold > 0:
        # disable highlighting for short rearrangements
        synteny_dict[genome].loc[(synteny_dict[genome][query_block_len_column_name] < args.min_len_threshold) & (synteny_dict[genome][target_block_len_column_name] < args.min_len_threshold),
                                 "connector_color"] = args.default_color
        synteny_dict[genome].loc[(synteny_dict[genome][query_block_len_column_name] < args.min_len_threshold) & (synteny_dict[genome][target_block_len_column_name] < args.min_len_threshold),
                                 "connector_zorder"] = 0

#--------------------------------------------------------------------------------------


for index, genome in zip(range(0, len(genome_orderlist) - 1), genome_orderlist[:-1]):
    synteny_dict[genome].to_csv("{0}.{1}.to.{2}.raw.tab".format(args.output_prefix,
                                                                genome,
                                                                genome_orderlist[index + 1]),
                                sep="\t", index=False, header=True)
#----------------------------- Renaming of scaffolds ----------------------------------
for genome_index in range(0, len(genome_orderlist)-1):
    synteny_dict[genome_orderlist[genome_index]][target_scaffold_id_column_name].replace(syn_df_dict[genome_orderlist[genome_index + 1]]["syn"].to_dict(),
                                                                                         inplace=True)
    synteny_dict[genome_orderlist[genome_index]][query_scaffold_id_column_name].replace(syn_df_dict[genome_orderlist[genome_index]]["syn"].to_dict(),
                                                                                        inplace=True)
    synteny_dict[genome_orderlist[genome_index]].sort_values(by=[query_scaffold_id_column_name,
                                                                 query_start_column_name,
                                                                 query_end_column_name,
                                                                 target_scaffold_id_column_name,
                                                                 target_start_column_name,
                                                                 target_end_column_name],
                                                             inplace=True)
#for index, genome in zip(range(0, len(genome_orderlist) - 1), genome_orderlist[:-1]):
#    synteny_dict[genome].to_csv("{0}.{1}.to.{2}.raw.renamed.tab".format(args.output_prefix,
#                                                                        genome,
#                                                                        genome_orderlist[index + 1]),
#                                sep="\t", index=False, header=True)

#--------------------------------------------------------------------------------------
#-------------------------- Inversion of coordinates ----------------------------------
for genome_index in range(0, genome_number):
    genome = genome_orderlist[genome_index]

    print("Inverting (if necessary) {0} scaffolds...".format(genome))
    #print(centromere_df_dict[genome])
    centromere_df_dict[genome] = invert_coordinates_in_region_table(centromere_df_dict[genome], invertlist_series_dict[genome],
                                                                    lenlist_df_dict[genome], "scaffold_id", "start", "end",
                                                                    inverted_scaffolds_label="'")

    if genome_index < (genome_number - 1):  # apply listed inversion for all query genomes, i.e. for all except the last genome
        print("Inverting query coordinates in synteny file...")
        #print(lenlist_df_dict[genome])
        synteny_dict[genome] = invert_coordinates_in_synteny_table(synteny_dict[genome],
                                                                   invertlist_series_dict[genome],
                                                                   lenlist_df_dict[genome],
                                                                   query_scaffold_id_column_name,
                                                                   query_start_column_name,
                                                                   query_end_column_name,
                                                                   strand_column_name,
                                                                   inverted_scaffold_label)
    if genome_index > 0:  # apply listed inversion for all target genomes, i.e. for all except the first genome
        #print(genome_orderlist[genome_index - 1])
        #print(synteny_dict[genome_orderlist[genome_index - 1]])
        print("Inverting target coordinates in synteny file...")
        synteny_dict[genome_orderlist[genome_index - 1]] = invert_coordinates_in_synteny_table(synteny_dict[genome_orderlist[genome_index - 1]],
                                                                                               invertlist_series_dict[genome],
                                                                                               lenlist_df_dict[genome],
                                                                                               target_scaffold_id_column_name,
                                                                                               target_start_column_name,
                                                                                               target_end_column_name,
                                                                                               strand_column_name,
                                                                                               inverted_scaffold_label)
    lenlist_df_dict[genome].rename(index=dict(zip(invertlist_series_dict[genome], [scaf + inverted_scaffold_label for scaf in invertlist_series_dict[genome]])),
                                   inplace=True)
    #chr_color_df_dict[species].rename(index=dict(zip(invert_list_dict[species], [scaf + label  for scaf in invert_list_dict[species]])),inplace=True)
    orderlist_series_dict[genome].replace(dict(zip(invertlist_series_dict[genome], [scaf + inverted_scaffold_label for scaf in invertlist_series_dict[genome]])),
                                          inplace=True)
    color_df_dict[genome].rename(index=dict(zip(invertlist_series_dict[genome], [scaf + inverted_scaffold_label for scaf in invertlist_series_dict[genome]])),
                                 inplace=True)
#for genome in genome_orderlist:




#--------------------------------------------------------------------------------------


inversion_hex_color = mpl.colors.cnames[args.inversion_color]
translocation_hex_color = mpl.colors.cnames[args.translocation_color]
default_hex_color = mpl.colors.cnames["lightgrey"]
long_block_hex_color = "#00FF00"
short_block_hex_color = "#F4EA56"

for index, genome in zip(range(0, len(genome_orderlist) - 1), genome_orderlist[:-1]):
    target_genome = genome_orderlist[index + 1]
    output_pr = "{0}.{1}.to.{2}".format(args.output_prefix, genome, target_genome)

    synteny_dict[genome].sort_values(by=[query_scaffold_id_column_name,
                                         query_start_column_name,
                                         query_end_column_name,
                                         target_scaffold_id_column_name,
                                         target_start_column_name,
                                         target_end_column_name], inplace=True)

    synteny_dict[genome].to_csv(output_pr + ".tab", sep="\t", index=False, header=True)
    writer = pd.ExcelWriter('{0}.xlsx'.format(output_pr), engine='xlsxwriter')
    workbook = writer.book
    # Adjust default format
    sheet_name = "{0}.to.{1}".format(genome, target_genome)

    green_hex = "#00FF00"
    light_green_hex = "#90EE90"
    light_blue_hex = "#ADD8E6"
    light_yellow_hex = "#F4EA56"
    light_orange_hex = "#FFD580"
    light_red_hex = "#FF7377"
    long_block_format = workbook.add_format({'bg_color': light_green_hex})
    short_block_format = workbook.add_format({'bg_color': light_blue_hex})
    too_short_block_format = workbook.add_format({'bg_color': light_orange_hex})

    inversion_format = workbook.add_format({'bg_color': inversion_hex_color})
    translocation_format = workbook.add_format({'bg_color': translocation_hex_color})
    default_format = workbook.add_format({'bg_color': default_hex_color})
    long_block_format = workbook.add_format({'bg_color': long_block_hex_color})
    short_block_format = workbook.add_format({'bg_color': short_block_hex_color})

    type_column_name = "type"
    connector_column_name = "connector_color"
    #print(bed_dict)
    #for species in bed_dict:
    #    species_format_dict[species] = {}
    #    for scaffold in color_df_dict[species].index:
    #        species_format_dict[species][scaffold] = workbook.add_format(
    #            {'bg_color': color_df_dict[species].loc[scaffold, "color"]})
    row_number = len(synteny_dict[genome])
    column_number = len(synteny_dict[genome].columns)
    query_column_idx = list(synteny_dict[genome].columns).index(query_scaffold_id_column_name)
    target_column_idx = list(synteny_dict[genome].columns).index(target_scaffold_id_column_name)
    query_block_len_column_idx = list(synteny_dict[genome].columns).index(query_block_len_column_name)
    target_block_len_column_idx = list(synteny_dict[genome].columns).index(target_block_len_column_name)
    type_column_idx = list(synteny_dict[genome].columns).index(type_column_name)
    connector_color_column_idx = list(synteny_dict[genome].columns).index(connector_column_name)

    #print(color_df_dict[genome])
    query_scaffold_format_dict = {}
    target_scaffold_format_dict = {}
    for scaffold_id in color_df_dict[genome].index:
        query_scaffold_format_dict[scaffold_id] = workbook.add_format({'bg_color': color_df_dict[genome].loc[scaffold_id, "color"]})

    for scaffold_id in color_df_dict[target_genome].index:
        target_scaffold_format_dict[scaffold_id] = workbook.add_format({'bg_color': color_df_dict[target_genome].loc[scaffold_id, "color"]})

    synteny_dict[genome].to_excel(writer, sheet_name=sheet_name,
                                  freeze_panes=(1, 2), index=False)

    for row in range(1, row_number + 1):
        if synteny_dict[genome][query_scaffold_id_column_name].iloc[row - 1] in query_scaffold_format_dict:
            writer.sheets[sheet_name].write(row, query_column_idx, synteny_dict[genome][query_scaffold_id_column_name].iloc[row - 1],
                                             # color query column
                                              query_scaffold_format_dict[synteny_dict[genome][query_scaffold_id_column_name].iloc[row - 1]])

        if synteny_dict[genome][target_scaffold_id_column_name].iloc[row - 1] in target_scaffold_format_dict:
            writer.sheets[sheet_name].write(row, target_column_idx, synteny_dict[genome][target_scaffold_id_column_name].iloc[row - 1],
                                             # color query column
                                            target_scaffold_format_dict[synteny_dict[genome][target_scaffold_id_column_name].iloc[row - 1]])

        writer.sheets[sheet_name].write(row, query_block_len_column_idx, synteny_dict[genome][query_block_len_column_name].iloc[row - 1],
                                        short_block_format if synteny_dict[genome][query_block_len_column_name].iloc[row - 1] < args.min_len_threshold else long_block_format)
        writer.sheets[sheet_name].write(row, target_block_len_column_idx, synteny_dict[genome][target_block_len_column_name].iloc[row - 1],
                                        short_block_format if synteny_dict[genome][target_block_len_column_name].iloc[row - 1] < args.min_len_threshold else long_block_format)

        writer.sheets[sheet_name].write(row, type_column_idx, synteny_dict[genome][type_column_name].iloc[row - 1],
                                        inversion_format if synteny_dict[genome][type_column_name].iloc[row - 1] == "inversion" else translocation_format if synteny_dict[genome][type_column_name].iloc[row - 1] == "translocation" else default_format)
        writer.sheets[sheet_name].write(row, connector_color_column_idx, synteny_dict[genome][connector_column_name].iloc[row - 1],
                                        inversion_format if synteny_dict[genome][connector_column_name].iloc[row - 1] == args.inversion_color else translocation_format if synteny_dict[genome][connector_column_name].iloc[row - 1] == args.translocation_color else default_format)
    workbook.worksheets()[0].autofit()
    workbook.formats[0].set_align('center')
    writer.close()



border_offset_fraction = 0.05
interchr_space_fraction = 0.3

maximal_x = max_genome_length * (1 + interchr_space_fraction)

height = 9
length = 1000000
distance = args.genome_distance

maximal_y = height * len(lenlist_df_dict) + distance * (len(lenlist_df_dict) - 1)

x_start = 0
y_start = 0

zorder_dict = {
               "background": 1,
               "connector": 500,
               "chromosomes": 1000,
               "chromosome_lines": 1500,
               "label": 2000
               }

x_scale_factor = 1

fig = plt.figure(1, figsize=(args.figure_width, genome_number * args.figure_height_per_genome), dpi=300)
ax = plt.subplot()

ax.set_axis_off()

xmin = -4*border_offset_fraction * maximal_x
xmax = (1 + border_offset_fraction) * maximal_x
ymin = -9*border_offset_fraction * maximal_y
ymax = (1 + 9 * border_offset_fraction) * maximal_y

plt.xlim(xmin=xmin, xmax=xmax)
plt.ylim(ymin=ymin, ymax=ymax)


ax.add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, color="white", alpha=1.0))


def patch_function(row): #centromere_df):
    #print(row)
    #print("AAAAA")
    #print(pd.isna(row.iloc[5]))
    #print("BBBBBBBBbb")
    #print((row.iloc[5] + row.iloc[1]) if not pd.isna(row.iloc[5]) else None)
    #print((row.iloc[6] + row.iloc[1]) if not pd.isna(row.iloc[5]) else None)
    return LinearChromosome(row.iloc[1], row.iloc[2], row.iloc[0], height, rounded=True,
                            x_scale_factor=args.smooth_multiplicator * maximal_x/10 / maximal_y, #maximal_x/10 / maximal_y,
                            zorder=zorder_dict["chromosomes"],
                            edgecolor=row.iloc[3],
                            facecolor=row.iloc[3],
                            alpha=0.9,
                            linewidth=0.3,
                            centromere_start=(row.iloc[5]) if not pd.isna(row.iloc[5]) else None, #+ row.iloc[1]
                            centromere_end=(row.iloc[6]) if not pd.isna(row.iloc[5]) else None, #+ row.iloc[1]
                            show_centromere=True)


def chromosome_line_function(row, height):
    return Line2D(xdata=(row.iloc[1], row.iloc[1] + row.iloc[0]),
                  ydata=(row.iloc[2] + height/2, row.iloc[2] + height/2,),
                  color="black",
                  zorder=zorder_dict["chromosome_lines"],
                  alpha=0.4,
                  linewidth=0.3,)


color_number = len(genome_orderlist)
colors = distinctipy.get_colors(color_number)
color_list = list(map(rgb_tuple_to_hex, colors))

for species, index, color, species_label in zip(genome_orderlist, range(0, len(genome_orderlist)), color_list, genome_labellist): #genome_orderlist):
    #print(centromere_df_dict[species])
    interchr_space = ((maximal_x - total_len_dict[species]) / (chr_number_dict[species] - 1)) if chr_number_dict[species] > 1 else 0
    lenlist_df_dict[species]["x_offset"] = lenlist_df_dict[species]["length"].cumsum().shift(periods=1, fill_value=0) + np.array(range(0, chr_number_dict[species])) * interchr_space
    lenlist_df_dict[species]["y_offset"] = (height + distance) * index
    lenlist_df_dict[species]["color"] = color

    lenlist_df_dict[species]["label"] = pd.Series(list(lenlist_df_dict[species].index),
                                                  index=lenlist_df_dict[species].index).apply(lambda s: s[args.scaffold_prefix_cut:])
    lenlist_df_dict[species]["centromere_start"] = centromere_df_dict[species]["start"]
    lenlist_df_dict[species]["centromere_end"] = centromere_df_dict[species]["end"]
    #print(lenlist_df_dict[species])
    patch_collection = PatchCollection(lenlist_df_dict[species].apply(patch_function,
                                                                      axis=1),
                                       match_original=True,
                                       antialiased=False,
                                       zorder=zorder_dict["chromosomes"])
    ax.add_collection(patch_collection)

    for line in lenlist_df_dict[species].apply(partial(chromosome_line_function, height=height), axis=1):
        ax.add_line(line)

    if not args.hide_chromosome_labels:
        for chromosome in lenlist_df_dict[species].index:
            ax.annotate(lenlist_df_dict[species].loc[chromosome, "label"],
                        xy=(lenlist_df_dict[species].loc[chromosome, "x_offset"] + lenlist_df_dict[species].loc[chromosome, "length"]/2,
                            lenlist_df_dict[species].loc[chromosome, "y_offset"] + height), xycoords='data',
                        fontsize=args.chromosome_label_fontsize,
                        xytext=(0, 0), textcoords='offset points', rotation=args.chromosome_label_angle,
                        ha="center", va="bottom",
                        color="black",
                        zorder=zorder_dict["label"])

    ax.annotate(species_label,
                xy=(- border_offset_fraction / 2 * maximal_x, (height + distance) * index + height/2),
                xycoords='data',
                fontsize=args.genome_label_fontsize,
                fontstyle="italic",
                xytext=(0, 0), textcoords='offset points', rotation=args.genome_label_angle,
                ha="right", va="center",
                color="black",
                zorder=zorder_dict["label"])


def connector_function(row, length_df_dict, top_species, bottom_species, default_color="lightgrey", connector_color_idx=8,
                      top_scaffold_idx=3, top_start_idx=4, top_end_idx=5, bottom_scaffold_idx=0, bottom_start_idx=1, bottom_end_idx=2, strand_idx=6):
    con_len = (connector_color_idx + 1) if connector_color_idx is not None else None
    y_chr_shift = height / 2
    #print("AAAAAAAAA")
    #print(top_species)
    #print(length_df_dict[top_species])
    #print("BBBBBBBBB")
    #print(bottom_species)
    #print(length_df_dict[bottom_species])
    return CubicBezierConnector(
                                 (row.iloc[top_start_idx] + length_df_dict[top_species].loc[row.iloc[top_scaffold_idx], "x_offset"], length_df_dict[top_species].loc[row.iloc[top_scaffold_idx], "y_offset"] + y_chr_shift),
                                 (row.iloc[top_end_idx] + length_df_dict[top_species].loc[row.iloc[top_scaffold_idx], "x_offset"], length_df_dict[top_species].loc[row.iloc[top_scaffold_idx], "y_offset"] + y_chr_shift),

                                ((row.iloc[bottom_start_idx] if row.iloc[strand_idx] == "+" else row.iloc[bottom_end_idx]) + length_df_dict[bottom_species].loc[row.iloc[bottom_scaffold_idx], "x_offset"],
                                 length_df_dict[bottom_species].loc[row.iloc[bottom_scaffold_idx], "y_offset"] + y_chr_shift),

                                ((row.iloc[bottom_end_idx] if row.iloc[strand_idx] == "+" else row.iloc[bottom_start_idx]) + length_df_dict[bottom_species].loc[row.iloc[bottom_scaffold_idx], "x_offset"],
                                 length_df_dict[bottom_species].loc[row.iloc[bottom_scaffold_idx], "y_offset"] + y_chr_shift),

                                 x_fraction_parameter=2,
                                 y_fraction_parameter=2,
                                 y_shift=distance,
                                 edgecolor=default_color if (con_len is None) or (len(row) < con_len) else row.iloc[connector_color_idx] if row.iloc[connector_color_idx] != "default" else default_color,
                                 facecolor=default_color if (con_len is None) or (len(row) < con_len) else row.iloc[connector_color_idx] if row.iloc[connector_color_idx] != "default" else default_color,
                                 alpha=0.5 if (con_len is None) or (len(row) < con_len) else 1.0 if row.iloc[connector_color_idx] != "default" else 0.5,
                                 fill= True,
                                 zorder=zorder_dict["connector"]
                                 )


connector_collection_dict = {}

for genome, genome_index in zip(genome_orderlist[:-1], range(0, len(genome_orderlist) - 1)):
    if "connector_zorder" in synteny_dict[genome]:
        synteny_dict[genome]["connector_zorder"] += zorder_dict["connector"]
        connector_collection_dict[genome] = {}
        for zorder in sorted(synteny_dict[genome]["connector_zorder"].unique()):
            #print(lenlist_df_dict[genome])
            connector_collection_dict[genome][zorder] = PatchCollection(synteny_dict[genome][synteny_dict[genome]["connector_zorder"] == zorder].apply(partial(connector_function,
                                                                                                                                                                       length_df_dict=lenlist_df_dict,
                                                                                                                                                                       top_species=genome_orderlist[genome_index+1],
                                                                                                                                                                       connector_color_idx=connector_color_idx,
                                                                                                                                                                       bottom_species=genome,
                                                                                                                                                                       top_scaffold_idx=target_scaffold_idx,
                                                                                                                                                                       top_start_idx=target_start_idx,
                                                                                                                                                                       top_end_idx=target_end_idx,
                                                                                                                                                                       bottom_scaffold_idx=query_scaffold_idx,
                                                                                                                                                                       bottom_start_idx=query_start_idx,
                                                                                                                                                                       bottom_end_idx=query_end_idx,
                                                                                                                                                                       strand_idx=strand_idx), axis=1),
                                                                            match_original=True,
                                                                            antialiased=False,
                                                                            zorder=zorder)
            ax.add_collection(connector_collection_dict[genome][zorder])
    else:
        connector_collection_dict[genome] = PatchCollection(synteny_dict[genome].apply(partial(connector_function,
                                                                                                       length_df_dict=lenlist_df_dict,
                                                                                                       top_species=genome_orderlist[genome_index+1],
                                                                                                       connector_color_idx=connector_color_idx,
                                                                                                       bottom_species=genome,
                                                                                                       top_scaffold_idx=target_scaffold_idx,
                                                                                                       top_start_idx=target_start_idx,
                                                                                                       top_end_idx=target_end_idx,
                                                                                                       bottom_scaffold_idx=query_scaffold_idx,
                                                                                                       bottom_start_idx=query_start_idx,
                                                                                                       bottom_end_idx=query_end_idx,
                                                                                                       strand_idx=strand_idx), axis=1),
                                                            match_original=True,
                                                            antialiased=False,
                                                            zorder=zorder_dict["connector"])
        ax.add_collection(connector_collection_dict[genome])

plt.subplots_adjust(left=args.subplots_adjust_left, right=args.subplots_adjust_right, bottom=args.subplots_adjust_bottom,
                    top=args.subplots_adjust_top)
plt.title(args.title, fontsize=args.title_fontsize)
for ext in args.output_formats:
    plt.savefig("{0}.{1}".format(args.output_prefix, ext))
