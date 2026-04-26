#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os

import argparse
import glob
from functools import partial
from pathlib import Path
from collections import OrderedDict

import pandas as pd
import numpy as np
from copy import deepcopy
from distinctipy import distinctipy

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

from MACE.Routines import Synteny, Parsing
from MACE.Routines import Parsing
from MACE.Routines import Visualization
from MACE.Visualization.Polygons import LinearChromosome
from MACE.Visualization.Connectors import CubicBezierConnector
from RouToolPa.Parsers.PSL import CollectionPSL


def split_comma_separated_list(string):
    return string.split(",")


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
parser.add_argument("--genome_config", action="store", dest="genome_config", default=None,
                    help="Optional configuration file. If set orderlists and invertlists from genome folders will be ignored")

parser.add_argument("--strand_switch_label", action="store", dest="strand_switch_label", default="*",
                    help="Symbol to be used in the genome config as a strand switch label. "
                         "If in scaffold id there is only one such a symbol a query strand switch will be applied to the scaffold. "
                         "If two symbol - target strand switch, if three - both query and target strand switches, respectively. "
                         "Default: '*', i.e use '*' - for query strand switch, '**' - for target strand switch, '***' - for both. "
                         "Note, strand switch symbols must occupy the very last positions in the scaffold id, but BEFORE the inversion symbol. "
                         "If you wish also to invert the scaffold, inversion symbol must be the very last symbol in the scaffold id, ")

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
parser.add_argument("--inverted_scaffold_label", action="store", dest="inverted_scaffold_label", default="'",
                    help="Symbol to use for labeling inverted scaffolds. Must be a very last symbol in the scaffold id. Default: '")
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
                    help="Comma-separated list of formats (supported by matplotlib) of "
                         "output figure.Default: svg,png")

parser.add_argument("-l", "--title", action="store", dest="title", default="Macrosynteny",
                    help="Suptitle of figure. Default: 'Macrosynteny'")
parser.add_argument("--title_fontsize", action="store", dest="title_fontsize", default=20, type=int,
                    help="Fontsize of the figure. Default: 20")
parser.add_argument("--chromosome_height", action="store", dest="chromosome_height", default=9, type=float,
                    help="Height of chromosomes on the plot. Increase or decrease this parameter to make chromosomes "
                         "thicker or thinner. Default: 9")
parser.add_argument("--hide_chromosome_labels", action="store_true", dest="hide_chromosome_labels", default=False,
                    help="Hide chromosome labels. Default: False")

parser.add_argument("--manual_figure_adjustment", action="store_true", dest="manual_figure_adjustment", default=False,
                    help="Adjust borders of figure manually using options below. Default: False, i.e. scaling is done automatically.")
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
parser.add_argument("--remove_scaffolds_absent_in_orderlist", action="store_true",
                    dest="remove_scaffolds_absent_in_orderlist",
                    default=False,
                    help="Remove scaffolds absent in orderlist. Default: False")
parser.add_argument("--do_not_highlight_translocations", action="store_true",
                    dest="do_not_highlight_translocations",
                    default=False,
                    help="Do not highlight translocations. Default: False, i.e highlight")
parser.add_argument("--do_not_highlight_inversions", action="store_true",
                    dest="do_not_highlight_inversions",
                    default=False,
                    help="Do not highlight inversions. Default: False, i.e highlight")

args = parser.parse_args()

data_dir = args.input_dir
data_dir_path = Path(data_dir)

genome_orderlist = args.genome_orderlist[::-1] if args.invert_genome_order else args.genome_orderlist
if args.genome_labellist is None:
    genome_labellist = genome_orderlist
else:
    genome_labellist = args.genome_labellist[::-1] if args.invert_genome_order else args.genome_labellist

if args.do_not_highlight_translocations:
    args.translocation_color = "default"
if args.do_not_highlight_inversions:
    args.inversion_color = "default"

genome_auxiliary_dict = {"genomes": OrderedDict(),
                         "genome_colors": OrderedDict()}

for genome in genome_orderlist:
    # nonempty whitelist file is necessary for each genome
    # nonempty lenlist file is necessary for each genome
    genome_auxiliary_dict["genomes"][genome] = \
                      Parsing.read_mace_auxiliary_input(len_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".len"]),
                                                        whitelist_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".whitelist"]),
                                                        orderlist_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".orderlist"]),
                                                        invertlist_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".invertlist"]),
                                                        inverted_scaffold_label=args.inverted_scaffold_label,
                                                        queryswitchstrandlist_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".queryswitchstrandlist"]),
                                                        targetswitchstrandlist_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".targetswitchstrandlist"]),
                                                        syn_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".syn"]),
                                                        syn_file_key_column=args.syn_file_key_column, syn_file_value_column=args.syn_file_value_column,
                                                        centromere_bed=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".centromere.bed"]),
                                                        scaffold_color_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".scaffold.color"]),
                                                        )

# Read genome config if it was set
Parsing.update_genome_dict_from_genome_config(genome_auxiliary_dict, args.genome_config,
                                              strand_switch_label=args.strand_switch_label,
                                              inverted_scaffold_label=args.inverted_scaffold_label)

# -------------------------------------------------------------------------------------------------------

#filter len list
for genome in genome_orderlist:     # genome_auxiliary_dict["genomes"][genome]["len_df"]
    genome_auxiliary_dict["genomes"][genome]["len_df"] = genome_auxiliary_dict["genomes"][genome]["len_df"].loc[genome_auxiliary_dict["genomes"][genome]["len_df"].index.isin(genome_auxiliary_dict["genomes"][genome]["whitelist_series"])]
    if genome_auxiliary_dict["genomes"][genome]["syn_dict"]:
        genome_auxiliary_dict["genomes"][genome]["len_df"].rename(index=genome_auxiliary_dict["genomes"][genome]["syn_dict"], inplace=True)

    genome_auxiliary_dict["genomes"][genome]["len_df"] = genome_auxiliary_dict["genomes"][genome]["len_df"].reindex(genome_auxiliary_dict["genomes"][genome]["orderlist_series"]).dropna()

# count number of chromosomes/scaffold and total length of them
total_len_dict = {genome: sum(genome_auxiliary_dict["genomes"][genome]["len_df"]["length"]) for genome in genome_orderlist}
chr_number_dict = {genome: len(genome_auxiliary_dict["genomes"][genome]["len_df"]) for genome in genome_orderlist}

max_genome_length = max(list(total_len_dict.values()))

if args.synteny_format == "psl":
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
        print("Filename: {0}".format(Parsing.get_filenames_for_extension(data_dir_path / genome_orderlist[genome_index],
                                                                         extension_list=["psl", "psl.gz"])))
        print("Query: {0}".format(genome_orderlist[genome_index]))
        print("Target: {0}".format(genome_orderlist[genome_index + 1]))
        synteny_dict[genome_orderlist[genome_index]] = CollectionPSL(Parsing.get_filenames_for_extension(data_dir_path / genome_orderlist[genome_index],
                                                                                                         extension_list=["psl", "psl.gz"]),
                                                                     target_white_list=genome_auxiliary_dict["genomes"][genome_orderlist[genome_index + 1]]["whitelist_series"],
                                                                     query_white_list=genome_auxiliary_dict["genomes"][genome_orderlist[genome_index]]["whitelist_series"],
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
#-------------------------------------------------------------------------------------------------
genome_number = len(genome_orderlist)

#------------------ Preprocessing of synteny blocks -----------------------------------
##----------------------- Detect nested blocks ----------------------------------------
detect_nested_blocks_preset = partial(Synteny.detect_nested_blocks,
                                      nested_in_block_column_name="nested_in",
                                      query_nested_in_block_column_name="query_nested_in",
                                      query_scaffold_id_column_name=query_scaffold_id_column_name,
                                      query_start_column_name=query_start_column_name,
                                      query_end_column_name=query_end_column_name,
                                      target_nested_in_block_column_name="target_nested_in",
                                      target_scaffold_id_column_name=target_scaffold_id_column_name,
                                      target_start_column_name=target_start_column_name,
                                      target_end_column_name=target_end_column_name)

tmp_dict = {}
block_remove_dict = {}
for genome_index in range(0, genome_number - 1):
    genome = genome_orderlist[genome_index]
    block_remove_dict[genome] = []
    tmp_dict[genome] = Synteny.detect_same_coords_blocks(synteny_dict[genome],
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

    if args.remove_same_coords_blocks and not tmp_dict[genome].empty:
        block_remove_dict[genome] = list(tmp_dict[genome]["target_same_coords"].dropna())

    if args.remove_nested_blocks and not tmp_dict[genome].empty:
        block_remove_dict[genome] += list(tmp_dict[genome][tmp_dict[genome]["query_nested_in"].notna()]["synteny_block_id"])
        block_remove_dict[genome] += list(tmp_dict[genome][tmp_dict[genome]["target_nested_in"].notna()]["synteny_block_id"])

    block_remove_dict[genome] = set(block_remove_dict[genome])
    synteny_dict[genome] = synteny_dict[genome][~synteny_dict[genome]["synteny_block_id"].isin(block_remove_dict[genome])]

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

    synteny_dict[genome]["type"] = "normal"
    synteny_dict[genome]["connector_color"] = args.default_color
    synteny_dict[genome]["connector_zorder"] = 0
    synteny_dict[genome].set_index([query_scaffold_id_column_name, target_scaffold_id_column_name], inplace=True)

    # major_strand_series is a series with two level index(qName, tName)
    major_strand_series = hit_sum.groupby(by=[query_scaffold_id_column_name,
                                              target_scaffold_id_column_name]).apply(lambda df: df.idxmax()[query_block_len_column_name][2])

    minus_strand_index = major_strand_series == "-"
    plus_strand_index = major_strand_series == "+"

    if args.invert_major_strand:  # first apply a global major strand switch

        print("Flag --invert_major_strand flag was set. Switching major strand for all genomes and all chromosomes...")
        major_strand_series[minus_strand_index] = "+"
        major_strand_series[plus_strand_index] = "-"

    # apply a local major strand switch
    if not genome_auxiliary_dict["genomes"][genome]["queryswitchstrandlist_series"].empty:
        print(f"Switching major strand for {genome} as query for {genome} vs {target_genome} alignment...")

        switch_series = pd.Series(major_strand_series.index.get_level_values(query_scaffold_id_column_name).isin(genome_auxiliary_dict["genomes"][genome]["queryswitchstrandlist_series"]))
        switch_series.index = major_strand_series.index

        minus_strand_index = switch_series & (major_strand_series == "-")
        plus_strand_index = switch_series & (major_strand_series == "+")
        major_strand_series[minus_strand_index] = "+"
        major_strand_series[plus_strand_index] = "-"

    # apply switchstrand list of target genome.
    if not genome_auxiliary_dict["genomes"][target_genome]["targetswitchstrandlist_series"].empty:
        print(f"Switching major strand for {target_genome} as target for {genome} vs {target_genome} alignment...")
        switch_series = pd.Series(major_strand_series.index.get_level_values(target_scaffold_id_column_name).isin(genome_auxiliary_dict["genomes"][target_genome]["targetswitchstrandlist_series"]))
        switch_series.index = major_strand_series.index

        minus_strand_index = switch_series & (major_strand_series == "-")
        plus_strand_index = switch_series & (major_strand_series == "+")
        major_strand_series[minus_strand_index] = "+"
        major_strand_series[plus_strand_index] = "-"

    synteny_dict[genome]["major_strand"] = major_strand_series
    synteny_dict[genome].reset_index(level=1, drop=False, inplace=True)
    major_query_homolog_series = hit_sum.droplevel(level=2).groupby(by=[query_scaffold_id_column_name,
                                                                        target_scaffold_id_column_name]).sum().groupby(by=[query_scaffold_id_column_name]).apply(lambda df: df.idxmax()[query_block_len_column_name][1])

    synteny_dict[genome]["major_query_homolog"] = major_query_homolog_series
    synteny_dict[genome].loc[synteny_dict[genome][target_scaffold_id_column_name] != synteny_dict[genome]["major_query_homolog"], "connector_color"] = args.translocation_color
    synteny_dict[genome].loc[synteny_dict[genome][target_scaffold_id_column_name] != synteny_dict[genome]["major_query_homolog"], "type"] = "translocation"
    synteny_dict[genome].loc[synteny_dict[genome][target_scaffold_id_column_name] != synteny_dict[genome]["major_query_homolog"], "connector_zorder"] = 50

    #detect translocations from target side
    synteny_dict[genome].reset_index(level=0, drop=False, inplace=True)
    major_target_homolog_series = hit_sum.droplevel(level=2).groupby(by=[target_scaffold_id_column_name,
                                                                         query_scaffold_id_column_name]).sum().groupby(by=[target_scaffold_id_column_name]).apply(lambda df: df.idxmax()[query_block_len_column_name][1])

    synteny_dict[genome].set_index(target_scaffold_id_column_name, inplace=True)
    synteny_dict[genome]["major_target_homolog"] = major_target_homolog_series
    synteny_dict[genome].loc[synteny_dict[genome][query_scaffold_id_column_name] != synteny_dict[genome]["major_target_homolog"], "connector_color"] = args.translocation_color
    synteny_dict[genome].loc[synteny_dict[genome][query_scaffold_id_column_name] != synteny_dict[genome]["major_target_homolog"], "type"] = "translocation"
    synteny_dict[genome].loc[synteny_dict[genome][query_scaffold_id_column_name] != synteny_dict[genome]["major_target_homolog"], "connector_zorder"] = 50

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
    synteny_dict[genome_orderlist[genome_index]][target_scaffold_id_column_name] = synteny_dict[genome_orderlist[genome_index]][target_scaffold_id_column_name].replace(genome_auxiliary_dict["genomes"][genome_orderlist[genome_index + 1]]["syn_dict"])
    synteny_dict[genome_orderlist[genome_index]][query_scaffold_id_column_name] = synteny_dict[genome_orderlist[genome_index]][query_scaffold_id_column_name].replace(genome_auxiliary_dict["genomes"][genome_orderlist[genome_index]]["syn_dict"])
    synteny_dict[genome_orderlist[genome_index]] = synteny_dict[genome_orderlist[genome_index]].sort_values(by=[query_scaffold_id_column_name,
                                                                                                            query_start_column_name,
                                                                                                            query_end_column_name,
                                                                                                            target_scaffold_id_column_name,
                                                                                                            target_start_column_name,
                                                                                                            target_end_column_name])
#--------------------------------------------------------------------------------------

if args.remove_scaffolds_absent_in_orderlist:
    for index in range(0, len(genome_orderlist) - 1):
        query_genome = genome_orderlist[index]
        target_genome = genome_orderlist[index + 1]
        bool_series = synteny_dict[query_genome][query_scaffold_id_column_name].isin(genome_auxiliary_dict["genomes"][query_genome]["orderlist_series"])
        bool_series &= synteny_dict[query_genome][target_scaffold_id_column_name].isin(genome_auxiliary_dict["genomes"][target_genome]["orderlist_series"])
        synteny_dict[query_genome] = synteny_dict[query_genome][bool_series]
#-------------------------- Inversion of coordinates ----------------------------------
for genome_index in range(0, genome_number):
    genome = genome_orderlist[genome_index]

    print("Inverting (if necessary) {0} scaffolds...".format(genome))
    genome_auxiliary_dict["genomes"][genome]["centromere_df"] = Parsing.invert_coordinates_in_region_table(genome_auxiliary_dict["genomes"][genome]["centromere_df"],
                                                                                                           genome_auxiliary_dict["genomes"][genome]["invertlist_series"],
                                                                                                           genome_auxiliary_dict["genomes"][genome]["len_df"],
                                                                                                           "scaffold", "start", "end",
                                                                                                           inverted_scaffolds_label=args.inverted_scaffold_label)

    if genome_index < (genome_number - 1):  # apply listed inversion for all query genomes, i.e. for all except the last genome
        print("Inverting query coordinates in synteny file...")
        synteny_dict[genome] = Parsing.invert_coordinates_in_synteny_table(synteny_dict[genome],
                                                                           genome_auxiliary_dict["genomes"][genome]["invertlist_series"],
                                                                           genome_auxiliary_dict["genomes"][genome]["len_df"],
                                                                           query_scaffold_id_column_name,
                                                                           query_start_column_name,
                                                                           query_end_column_name,
                                                                           strand_column_name,
                                                                           args.inverted_scaffold_label)
    if genome_index > 0:  # apply listed inversion for all target genomes, i.e. for all except the first genome
        print("Inverting target coordinates in synteny file...")
        synteny_dict[genome_orderlist[genome_index - 1]] = Parsing.invert_coordinates_in_synteny_table(synteny_dict[genome_orderlist[genome_index - 1]],
                                                                                                       genome_auxiliary_dict["genomes"][genome]["invertlist_series"],
                                                                                                       genome_auxiliary_dict["genomes"][genome]["len_df"],
                                                                                                       target_scaffold_id_column_name,
                                                                                                       target_start_column_name,
                                                                                                       target_end_column_name,
                                                                                                       strand_column_name,
                                                                                                       args.inverted_scaffold_label)
    genome_auxiliary_dict["genomes"][genome]["len_df"].rename(index=dict(zip(genome_auxiliary_dict["genomes"][genome]["invertlist_series"],
                                                                             [scaf + args.inverted_scaffold_label for scaf in genome_auxiliary_dict["genomes"][genome]["invertlist_series"]])),
                                                              inplace=True)
    genome_auxiliary_dict["genomes"][genome]["orderlist_series"].replace(dict(zip(genome_auxiliary_dict["genomes"][genome]["invertlist_series"],
                                                                                  [scaf + args.inverted_scaffold_label for scaf in genome_auxiliary_dict["genomes"][genome]["invertlist_series"]])),
                                                                         inplace=True)
    if genome_auxiliary_dict["genomes"][genome]["scaffold_color_df"] is not None:
        genome_auxiliary_dict["genomes"][genome]["scaffold_color_df"].rename(index=dict(zip(genome_auxiliary_dict["genomes"][genome]["invertlist_series"],
                                                                                        [scaf + args.inverted_scaffold_label for scaf in genome_auxiliary_dict["genomes"][genome]["invertlist_series"]])),
                                                                             inplace=True)

#--------------------------------------------------------------------------------------

default_hex_color = mpl.colors.cnames["lightgrey"]
inversion_hex_color = mpl.colors.cnames[args.inversion_color] if args.inversion_color != "default" else default_hex_color
translocation_hex_color = mpl.colors.cnames[args.translocation_color] if args.translocation_color != "default" else default_hex_color

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
    sheet_name = "{0}.to.{1}".format(genome, target_genome).replace("'", "") # ' is not allowed in the sheetname
    if len(sheet_name) > 31:
        print(f"WARNING!!! Excel has a hardlimit of 31 char for sheetname. '{sheet_name}' is longer. Cutting to 31 chars...")
        sheet_name = sheet_name[:31]

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

    row_number = len(synteny_dict[genome])
    column_number = len(synteny_dict[genome].columns)
    query_column_idx = list(synteny_dict[genome].columns).index(query_scaffold_id_column_name)
    target_column_idx = list(synteny_dict[genome].columns).index(target_scaffold_id_column_name)
    query_block_len_column_idx = list(synteny_dict[genome].columns).index(query_block_len_column_name)
    target_block_len_column_idx = list(synteny_dict[genome].columns).index(target_block_len_column_name)
    type_column_idx = list(synteny_dict[genome].columns).index(type_column_name)
    connector_color_column_idx = list(synteny_dict[genome].columns).index(connector_column_name)

    query_scaffold_format_dict = {}
    target_scaffold_format_dict = {}
    for scaffold_id in genome_auxiliary_dict["genomes"][genome]["scaffold_color_df"].index:
        query_scaffold_format_dict[scaffold_id] = workbook.add_format({'bg_color': genome_auxiliary_dict["genomes"][genome]["scaffold_color_df"].loc[scaffold_id, "color"]})

    for scaffold_id in genome_auxiliary_dict["genomes"][target_genome]["scaffold_color_df"].index:
        target_scaffold_format_dict[scaffold_id] = workbook.add_format({'bg_color': genome_auxiliary_dict["genomes"][target_genome]["scaffold_color_df"].loc[scaffold_id, "color"]})

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
    workbook.get_worksheet_by_name(sheet_name).autofit()
    workbook.formats[0].set_align('center')
    writer.close()

border_offset_fraction = 0.05
interchr_space_fraction = 0.3

maximal_x = max_genome_length * (1 + interchr_space_fraction)

height = args.chromosome_height
length = 1000000
distance = args.genome_distance

maximal_y = height * len(genome_auxiliary_dict["genomes"]) + distance * (len(genome_auxiliary_dict["genomes"]) - 1)

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


def patch_function(row):
    return LinearChromosome(row.iloc[1], row.iloc[2], row.iloc[0], height, rounded=True, #stranded=True,
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


colors = list(map(Visualization.rgb_tuple_to_hex,
                  distinctipy.get_colors(genome_number)))

if "genome_color_df" in genome_auxiliary_dict:
    for genome_index in range(0, genome_number):
        if pd.isna(genome_auxiliary_dict["genome_color_df"].loc[genome_orderlist[genome_index], "color"]):
            genome_auxiliary_dict["genome_color_df"].loc[genome, "color"] = colors[genome_index]
else:
    genome_auxiliary_dict["genome_color_df"] = pd.DataFrame.from_records(zip(genome_orderlist, colors), columns=["genome", "color"], index="genome")

genome_auxiliary_dict["genome_color_df"].to_csv(f"{args.output_prefix}.genome.color", sep="\t", index=False, header=True)

lenlist_df_dict = {genome: genome_auxiliary_dict["genomes"][genome]["len_df"] for genome in genome_auxiliary_dict["genomes"]}

for species, index, species_label in zip(genome_orderlist, range(0, len(genome_orderlist)), genome_labellist): #genome_orderlist):
    interchr_space = ((maximal_x - total_len_dict[species]) / (chr_number_dict[species] - 1)) if chr_number_dict[species] > 1 else 0
    lenlist_df_dict[species]["x_offset"] = lenlist_df_dict[species]["length"].cumsum().shift(periods=1, fill_value=0) + np.array(range(0, chr_number_dict[species])) * interchr_space
    lenlist_df_dict[species]["y_offset"] = (height + distance) * index
    lenlist_df_dict[species]["color"] = genome_auxiliary_dict["genome_color_df"].loc[species, "color"]

    lenlist_df_dict[species]["label"] = pd.Series(list(lenlist_df_dict[species].index),
                                                  index=lenlist_df_dict[species].index).apply(lambda s: s[args.scaffold_prefix_cut:])
    lenlist_df_dict[species]["centromere_start"] = genome_auxiliary_dict["genomes"][species]["centromere_df"]["start"]
    lenlist_df_dict[species]["centromere_end"] = genome_auxiliary_dict["genomes"][species]["centromere_df"]["end"]
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
                                 fill=True,
                                 zorder=zorder_dict["connector"]
                                 )


connector_collection_dict = {}

for genome, genome_index in zip(genome_orderlist[:-1], range(0, len(genome_orderlist) - 1)):
    if "connector_zorder" in synteny_dict[genome]:
        synteny_dict[genome]["connector_zorder"] += zorder_dict["connector"]
        connector_collection_dict[genome] = {}
        for zorder in sorted(synteny_dict[genome]["connector_zorder"].unique()):
            connector_collection_dict[genome][zorder] = \
                PatchCollection(synteny_dict[genome][synteny_dict[genome]["connector_zorder"] == zorder].apply(partial(connector_function,
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

plt.title(args.title, fontsize=args.title_fontsize)
if args.manual_figure_adjustment:
    plt.subplots_adjust(left=args.subplots_adjust_left, right=args.subplots_adjust_right, bottom=args.subplots_adjust_bottom,
                        top=args.subplots_adjust_top)

for ext in args.output_formats:
    if args.manual_figure_adjustment:
        plt.savefig(f"{args.output_prefix}.{ext}")
    else:
        plt.savefig(f"{args.output_prefix}.{ext}", bbox_inches="tight",)
