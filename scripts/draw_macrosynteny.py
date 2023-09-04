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
    temp_df.to_csv("tmp", sep="\t", index=True, header=True)
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


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True,
                    help="Input directory with data. Must contain subfolder for each genome. "
                         "Subfolders should have same name as genomes in --genome_orderlist"
                         "Each subfolder should contain: *.whitelist, *.len and synteny file "
                         "(except for the last genome). *.orderlist, *.invertlist and *.syn file are optional")
parser.add_argument("--genome_orderlist", action="store", dest="genome_orderlist", required=True,
                    type=split_comma_separated_list,
                    help="Comma-separated list of genomes to be used in figure.")
parser.add_argument("--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in synonym file. Default: 0")
parser.add_argument("--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in synonym file synonym. Default: 1")
parser.add_argument("--synteny_format", action="store", dest="synteny_format", default="psl",
                    help="Format of the synteny file. Allowed: psl(default)")

parser.add_argument("--min_len_threshold", action="store", dest="min_len_threshold", default=0, type=int,
                    help="Minimum length of rearranged block to be highlighted. "
                         "Default: 0, i.e. all rearranged blocks will be highlighted")

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

genome_orderlist = args.genome_orderlist
syn_file_key_column, syn_file_value_column = args.syn_file_key_column, args.syn_file_value_column

synteny_format = args.synteny_format

inverted_scaffold_label = "'"

# read files
print(data_dir_path)
# whitelist is obligatory
whitelist_series_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["whitelist"]),
                                             sep="\t", header=None, comment="#").squeeze("columns") for genome in genome_orderlist}
# orderlist might be absent in folders
orderlist_series_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["orderlist"]),
                                             sep="\t", header=None, comment="#").squeeze("columns") if get_filenames_for_extension(data_dir_path / genome, extension_list=["orderlist"]) is not None else pd.Series(dtype=str)  for genome in genome_orderlist}
# invertlist might be absent in folders
invertlist_series_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["invertlist"]),
                                             sep="\t", header=None, comment="#").squeeze("columns")  if get_filenames_for_extension(data_dir_path / genome, extension_list=["invertlist"]) is not None else pd.Series(dtype=str) for genome in genome_orderlist}
# lenlist is obligatory
lenlist_df_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["len"]),
                                             sep="\t", header=None, comment="#", names=["scaffold", "length"], index_col=0) for genome in genome_orderlist}
# synfile might be absent
syn_df_dict = {genome: pd.read_csv(get_filenames_for_extension(data_dir_path / genome, extension_list=["syn"]), usecols=(syn_file_key_column, syn_file_value_column),
                                             sep="\t", header=None, comment="#",
                                   names=["key", "syn"] if syn_file_key_column <= syn_file_value_column else ["syn", "key"]).set_index("key") if get_filenames_for_extension(data_dir_path / genome, extension_list=["syn"]) is not None else pd.DataFrame(columns=["syn"]) for genome in genome_orderlist}
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
    target_block_len_column_name = "qHitLen"

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
                                                                     target_syn_dict=syn_df_dict[genome_orderlist[genome_index + 1]]["syn"].to_dict() ,
                                                                     query_syn_dict=syn_df_dict[genome_orderlist[genome_index]]["syn"].to_dict(),
                                                                     parsing_mode="coordinates_only")
#print(syn_df_dict)
#print(synteny_dict["dingo"].records)
print("\n\n")
genome_number = len(genome_orderlist)
for genome_index in range(0, genome_number):
    genome = genome_orderlist[genome_index]
    print("Inverting (if necessary) {0} scaffolds...".format(genome))
    if genome_index < (genome_number - 1):  # apply listed inversion for all query genomes, i.e for all except the last genome
        print("Inverting query coordinates in synteny file...")
        #print(lenlist_df_dict[genome])
        synteny_dict[genome].records = invert_coordinates_in_synteny_table(synteny_dict[genome].records,
                                                                           invertlist_series_dict[genome],
                                                                           lenlist_df_dict[genome],
                                                                           query_scaffold_id_column_name,
                                                                           query_start_column_name,
                                                                           query_end_column_name,
                                                                           strand_column_name,
                                                                           inverted_scaffold_label)
    if genome_index > 0:  # apply listed inversion for all target genomes, i.e for all except the first genome
        #print(genome_orderlist[genome_index - 1])
        #print(synteny_dict[genome_orderlist[genome_index - 1]].records)
        print("Inverting target coordinates in synteny file...")
        synteny_dict[genome_orderlist[genome_index - 1]].records = invert_coordinates_in_synteny_table(synteny_dict[genome_orderlist[genome_index - 1]].records,
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

if synteny_format == "psl":

    for genome in genome_orderlist[:-1]: # all query genomes
        pd.set_option('display.max_rows', 40)
        columns_list = list(synteny_dict[genome].records.columns)
        hit_sum = synteny_dict[genome].records[[strand_column_name,
                                                query_scaffold_id_column_name,
                                                target_scaffold_id_column_name,
                                                query_block_len_column_name]].groupby(by=[query_scaffold_id_column_name,
                                                                                          target_scaffold_id_column_name,
                                                                                          strand_column_name]).sum()

        #print(synteny_dict[genome].records)
        synteny_dict[genome].records["connector_color"] = "default"
        synteny_dict[genome].records["connector_zorder"] = 0
        synteny_dict[genome].records.set_index(["qName", "tName"], inplace=True)
        #print(synteny_dict[genome].records)
        major_strand_series = hit_sum.groupby(by=[query_scaffold_id_column_name,
                                                  target_scaffold_id_column_name]).apply(lambda df: df.idxmax()[query_block_len_column_name][2])
        synteny_dict[genome].records["major_strand"] = major_strand_series
        synteny_dict[genome].records.reset_index(level=1, drop=False, inplace=True)
        #print(genome)
        #print(hit_sum)
        #detect translocations from query side
        major_query_homolog_series = hit_sum.droplevel(level=2).groupby(by=[query_scaffold_id_column_name,
                                                                            target_scaffold_id_column_name]).sum().groupby(by=[query_scaffold_id_column_name]).apply(lambda df: df.idxmax()[query_block_len_column_name][1])

        synteny_dict[genome].records["major_query_homolog"] = major_query_homolog_series
        synteny_dict[genome].records.loc[synteny_dict[genome].records[target_scaffold_id_column_name] != synteny_dict[genome].records["major_query_homolog"], "connector_color"] = "blue"
        synteny_dict[genome].records.loc[synteny_dict[genome].records[target_scaffold_id_column_name] != synteny_dict[genome].records["major_query_homolog"], "connector_zorder"] = 50

        #detect translocations from target side
        synteny_dict[genome].records.reset_index(level=0, drop=False, inplace=True)
        major_target_homolog_series = hit_sum.droplevel(level=2).groupby(by=[target_scaffold_id_column_name,
                                                                             query_scaffold_id_column_name]).sum().groupby(by=[target_scaffold_id_column_name]).apply(lambda df: df.idxmax()[target_block_len_column_name][1])
        #print(major_target_homolog_series)

        synteny_dict[genome].records.set_index(target_scaffold_id_column_name, inplace=True)
        synteny_dict[genome].records["major_target_homolog"] = major_target_homolog_series
        synteny_dict[genome].records.loc[synteny_dict[genome].records[query_scaffold_id_column_name] != synteny_dict[genome].records["major_target_homolog"], "connector_color"] = "blue"
        synteny_dict[genome].records.loc[synteny_dict[genome].records[query_scaffold_id_column_name] != synteny_dict[genome].records["major_target_homolog"], "connector_zorder"] = 50
        #synteny_dict[genome].records.to_csv("AAAAAAAAAAAAAA.{0}.tmp".format(genome), sep="\t")
        #print(synteny_dict[genome].records)

        synteny_dict[genome].records.reset_index(level=0, drop=False, inplace=True)

        synteny_dict[genome].records.loc[synteny_dict[genome].records[strand_column_name] != synteny_dict[genome].records["major_strand"], "connector_color"] = "red"
        synteny_dict[genome].records.loc[synteny_dict[genome].records[strand_column_name] != synteny_dict[genome].records["major_strand"], "connector_zorder"] = 20

        connector_color_idx = len(columns_list)
        connector_zorder_idx = len(columns_list) + 1
        synteny_dict[genome].records = synteny_dict[genome].records[columns_list + ["connector_color", "connector_zorder"]]

        if args.min_len_threshold > 0:
            # disable highlighting for short rearrangements
            synteny_dict[genome].records.loc[(synteny_dict[genome].records[query_block_len_column_name] < args.min_len_threshold) & (synteny_dict[genome].records[target_block_len_column_name] < args.min_len_threshold),
                                             "connector_color"] = "default"
            synteny_dict[genome].records.loc[(synteny_dict[genome].records[query_block_len_column_name] < args.min_len_threshold) & (synteny_dict[genome].records[target_block_len_column_name] < args.min_len_threshold),
                                             "connector_zorder"] = 0

for index, genome in zip(range(0, len(genome_orderlist) - 1), genome_orderlist[:-1]):
    synteny_dict[genome].records.sort_values(by=["qName", "qStart", "qEnd", "tName", "tStart", "tEnd"]).to_csv("{0}.{1}.to.{2}.tab".format(args.output_prefix,
                                                                                                                                           genome,
                                                                                                                                           genome_orderlist[index + 1]),
                                             sep="\t", index=False, header=True)

border_offset_fraction = 0.05
interchr_space_fraction = 0.3

maximal_x = max_genome_length * (1 + interchr_space_fraction)

height = 9
length = 1000000
distance = 100

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


def patch_function(row,):
    return LinearChromosome(row[1], row[2], row[0], height, rounded=True,
                            x_scale_factor=maximal_x/10 / maximal_y,
                            zorder=zorder_dict["chromosomes"],
                            edgecolor=row[3],
                            facecolor=row[3],
                            alpha=0.9,
                            linewidth=0.3)


def chromosome_line_function(row, height):
    return Line2D(xdata=(row[1], row[1] + row[0]),
                  ydata=(row[2] + height/2, row[2] + height/2,),
                  color="black",
                  zorder=zorder_dict["chromosome_lines"],
                  alpha=0.4,
                  linewidth=0.3,)


color_number = len(genome_orderlist)
colors = distinctipy.get_colors(color_number)
color_list = list(map(rgb_tuple_to_hex, colors))

for species, index, color, species_label in zip(genome_orderlist, range(0, len(genome_orderlist)), color_list, genome_orderlist):
    #print(species)
    interchr_space = ((maximal_x - total_len_dict[species]) / (chr_number_dict[species] - 1)) if chr_number_dict[species] > 1 else 0
    lenlist_df_dict[species]["x_offset"] = lenlist_df_dict[species]["length"].cumsum().shift(periods=1, fill_value=0) + np.array(range(0, chr_number_dict[species])) * interchr_space
    lenlist_df_dict[species]["y_offset"] = (height + distance) * index
    lenlist_df_dict[species]["color"] = color

    lenlist_df_dict[species]["label"] = pd.Series(list(lenlist_df_dict[species].index),
                                                  index=lenlist_df_dict[species].index).apply(lambda s:s[3:])

    patch_collection = PatchCollection(lenlist_df_dict[species].apply(patch_function, axis=1),
                                       match_original=True,
                                       antialiased=False,
                                       zorder=zorder_dict["chromosomes"])
    ax.add_collection(patch_collection)

    for line in lenlist_df_dict[species].apply(partial(chromosome_line_function, height=height), axis=1):
        ax.add_line(line)

    for chromosome in lenlist_df_dict[species].index:
        ax.annotate(lenlist_df_dict[species].loc[chromosome, "label"],
                    xy=(lenlist_df_dict[species].loc[chromosome, "x_offset"] + lenlist_df_dict[species].loc[chromosome, "length"]/2,
                        lenlist_df_dict[species].loc[chromosome, "y_offset"] + height), xycoords='data',
                    fontsize=5,
                    xytext=(0, 0), textcoords='offset points',
                    ha="center", va="bottom",
                    color="black",
                    zorder=zorder_dict["label"])

    ax.annotate(species_label,
                xy=(- border_offset_fraction / 2 * maximal_x, (height + distance) * index + height/2),
                xycoords='data',
                fontsize=9,
                fontstyle= "italic",
                xytext=(0, 0), textcoords='offset points',
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
                                 (row[top_start_idx] + length_df_dict[top_species].loc[row[top_scaffold_idx], "x_offset"], length_df_dict[top_species].loc[row[top_scaffold_idx], "y_offset"] + y_chr_shift),
                                 (row[top_end_idx] + length_df_dict[top_species].loc[row[top_scaffold_idx], "x_offset"], length_df_dict[top_species].loc[row[top_scaffold_idx], "y_offset"] + y_chr_shift),

                                ((row[bottom_start_idx] if row[strand_idx] == "+" else row[bottom_end_idx]) + length_df_dict[bottom_species].loc[row[bottom_scaffold_idx], "x_offset"],
                                 length_df_dict[bottom_species].loc[row[bottom_scaffold_idx], "y_offset"] + y_chr_shift),

                                ((row[bottom_end_idx] if row[strand_idx] == "+" else row[bottom_start_idx]) + length_df_dict[bottom_species].loc[row[bottom_scaffold_idx], "x_offset"],
                                 length_df_dict[bottom_species].loc[row[bottom_scaffold_idx], "y_offset"] + y_chr_shift),

                                 x_fraction_parameter=2,
                                 y_fraction_parameter=2,
                                 y_shift=distance,
                                 edgecolor=default_color if (con_len is None) or (len(row) < con_len) else row[connector_color_idx] if row[connector_color_idx] != "default" else default_color,
                                 facecolor=default_color if (con_len is None) or (len(row) < con_len) else row[connector_color_idx] if row[connector_color_idx] != "default" else default_color,
                                 alpha=0.5 if (con_len is None) or (len(row) < con_len) else 1.0 if row[connector_color_idx] != "default" else 0.5,
                                 fill= True,
                                 zorder=zorder_dict["connector"]
                                 )


connector_collection_dict = {}

for genome, genome_index in zip(genome_orderlist[:-1], range(0, len(genome_orderlist) - 1)):
    if "connector_zorder" in synteny_dict[genome].records:
        synteny_dict[genome].records["connector_zorder"] += zorder_dict["connector"]
        connector_collection_dict[genome] = {}
        for zorder in sorted(synteny_dict[genome].records["connector_zorder"].unique()):
            #print(lenlist_df_dict[genome])
            connector_collection_dict[genome][zorder] = PatchCollection(synteny_dict[genome].records[synteny_dict[genome].records["connector_zorder"] == zorder].apply(partial(connector_function,
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
        connector_collection_dict[genome] = PatchCollection(synteny_dict[genome].records.apply(partial(connector_function,
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
