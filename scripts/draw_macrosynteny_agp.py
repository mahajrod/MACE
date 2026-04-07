#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import glob
import yaml
import logging
import logging.config
import argparse

from functools import partial
from pathlib import Path
from collections import OrderedDict
from collections.abc import Mapping

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
from RouToolPa.Parsers.AGP import CollectionAGP


def split_comma_separated_list(string):
    return string.split(",")


def rgb_tuple_to_hex(rgb_tuple):
    if isinstance(rgb_tuple, str):
        return rgb_tuple
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
    #print(filelist)
    if force_uniq:
        if len(filelist) > 1:
            raise ValueError("Found more than one file with extensions: {0} in directory {1}".format(",".join(extension_list, str(dir_path))))
        else:
            return filelist[0]

    return filelist


def copy_absent_entries(input_dictionary, output_dictionary):
    for entry in input_dictionary:
        if entry not in output_dictionary:
            output_dictionary[entry] = deepcopy(input_dictionary[entry])
        else:
            if not isinstance(output_dictionary[entry], Mapping): # check if existing entry is not dictionary or dictionary like
                continue # exit from recursion
            copy_absent_entries(input_dictionary[entry], output_dictionary[entry])


def invert_coordinates_in_synteny_table(df, scaffold_list, length_df, scaffold_column, start_column, end_column, strand_column, inverted_scaffolds_label="'"):
    temp_df = deepcopy(df)
    columns_list = list(temp_df.columns)
    temp_df["length_column"] = 0
    temp_df.set_index(scaffold_column, inplace=True)
    #print (temp_df)
    #temp_df.to_csv("tmp", sep="\t", index=True, header=True)
    #print(scaffold_list)
    #print(length_df)
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
    #print(length_df)
    for scaffold in temp_df.index.unique():
        #print(scaffold)
        #print(temp_df)
        #print(length_df)
        #print(length_df)
        if scaffold in length_df.index:
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
                         "chromosome, please, use genome=specific switchstrandlist file.  Default: False")
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
                    help="Comma-separated list of formats (supported by matlotlib) of "
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
parser.add_argument("--scaffold_label_fontsize", action="store", dest="scaffold_label_fontsize", type=int, default=None,
                    help="Fontsize of scaffold labels. Default: 7")
parser.add_argument("--genome_label_fontsize", action="store", dest="genome_label_fontsize", type=int, default=None,
                    help="Fontsize of genome labels. Default: 11")
parser.add_argument("--scaffold_label_angle", action="store", dest="scaffold_label_angle", type=int, default=None,
                    help="Angle to rotate scaffold labels on the plot. Default: 0")
parser.add_argument("--genome_label_angle", action="store", dest="genome_label_angle", type=int, default=None,
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

parser.add_argument("--interscaf_space_fraction", action="store", default=0.3,
                    dest="interscaf_space_fraction",
                    help="Fraction of maximal genome length to use for spacers between scaffolds. Default: 0.3")

parser.add_argument("--log_level", action="store", default="INFO",
                    dest="log_level",
                    help="Logging level for stdout output. File {output_prefix}.full.log captures all log messages. "
                         "Allowed: INFO(default), DBG_SCR, DEBUG. Note that ERROR, WARNING and CRITICAL messages"
                         " are always written to stderr and files only.")

parser.add_argument("--settings_level", action="store", default="scaffold",
                    dest="settings_level",
                    help="Force specific settings level to be used. Allowed: 'scaffold' (default), 'genome', 'default'."
                         "'default' - top level settings override settings of individual genomes;"
                         "'genome'  - genome level settings override settings of individual scaffolds/superscaffolds, absent settings are filled by defaults;"
                         "'scaffold' - scaffold level settings are used, absent settings are filled from genome settings (if set) or defaults (if not set)")


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

# -------- Config loggers --------
SPACE_PER_TAB = 4
TAB = " " * SPACE_PER_TAB
TAB_NUM_PER_LVL = 1

LOW_LVL_DEBUG_SCRIPT_NUM = 12
LOW_LVL_DEBUG_SCRIPT_NAME = "LOW_DBG_SCR"
LOW_LVL_DEBUG_SCRIPT_METHOD_NAME = LOW_LVL_DEBUG_SCRIPT_NAME.lower()

DEBUG_SCRIPT_NUM = 15
DEBUG_SCRIPT_NAME = "DBG_SCR"
DEBUG_SCRIPT_METHOD_NAME = DEBUG_SCRIPT_NAME.lower()

WARNING_SCRIPT_NUM = 25
WARNING_SCRIPT_NAME = "WARN_SCR"
WARNING_SCRIPT_METHOD_NAME = WARNING_SCRIPT_NAME .lower()

# Adding custom level to capture debug messages related to the script only

def logForLevelDBG_SCR(self, message, *args, **kwargs):
    if self.isEnabledFor(DEBUG_SCRIPT_NUM):
        self._log(DEBUG_SCRIPT_NUM, message, args, **kwargs)


def logToRootDBG_SCR(message, *args, **kwargs):
    logging.log(DEBUG_SCRIPT_NUM, message, *args, **kwargs)


def logForLevelLOW_LVL_DBG_SCR(self, message, *args, **kwargs):
    if self.isEnabledFor(LOW_LVL_DEBUG_SCRIPT_NUM):
        self._log(LOW_LVL_DEBUG_SCRIPT_NUM, message, args, **kwargs)


def logToRootLOW_LVL_DBG_SCR(message, *args, **kwargs):
    logging.log(LOW_LVL_DEBUG_SCRIPT_NUM, message, *args, **kwargs)


def logForLevelWARNING_SCR(self, message, *args, **kwargs):
    if self.isEnabledFor(WARNING_SCRIPT_NUM):
        self._log(WARNING_SCRIPT_NUM, message, args, **kwargs)


def logToRootWARNING_SCR(message, *args, **kwargs):
    logging.log(WARNING_SCRIPT_NUM, message, *args, **kwargs)


logging.addLevelName(LOW_LVL_DEBUG_SCRIPT_NUM, LOW_LVL_DEBUG_SCRIPT_NAME)
setattr(logging, LOW_LVL_DEBUG_SCRIPT_NAME, LOW_LVL_DEBUG_SCRIPT_NUM)
setattr(logging.getLoggerClass(), LOW_LVL_DEBUG_SCRIPT_METHOD_NAME, logForLevelLOW_LVL_DBG_SCR)
setattr(logging, LOW_LVL_DEBUG_SCRIPT_METHOD_NAME, logToRootLOW_LVL_DBG_SCR)

logging.addLevelName(DEBUG_SCRIPT_NUM, DEBUG_SCRIPT_NAME)
setattr(logging, DEBUG_SCRIPT_NAME, DEBUG_SCRIPT_NUM)
setattr(logging.getLoggerClass(), DEBUG_SCRIPT_METHOD_NAME, logForLevelDBG_SCR)
setattr(logging, DEBUG_SCRIPT_METHOD_NAME, logToRootDBG_SCR)

logging.addLevelName(WARNING_SCRIPT_NUM, WARNING_SCRIPT_NAME)
setattr(logging, WARNING_SCRIPT_NAME, WARNING_SCRIPT_NUM)
setattr(logging.getLoggerClass(), WARNING_SCRIPT_METHOD_NAME, logForLevelWARNING_SCR)
setattr(logging, WARNING_SCRIPT_METHOD_NAME, logToRootWARNING_SCR)

log_config_dict = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
                   "main": {
                            "format": "%(asctime)s:%(levelname)-8s %(message)s", # "format": "%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
                            "datefmt": "%Y-%b-%d-%H:%M:%S",
                           },
                   "file": {
                            "format": "%(asctime)s:%(levelname)-14s %(message)s", # "format": "%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
                            "datefmt": "%Y-%b-%d-%H:%M:%S",
                           }
                  },
    "filters": {
        "warnings_and_below": {
            "()" : "__main__.filter_maker",
            "level": "WARNING"
        }
    },
    "handlers": {
        "stdout": {
            "class": "logging.StreamHandler",
            "level": "INFO",
            "formatter": "main",
            "stream": "ext://sys.stdout",
            "filters": ["warnings_and_below"]
        },
        "stderr": {
            "class": "logging.StreamHandler",
            "level": "WARNING",
            "formatter": "main",
            "stream": "ext://sys.stderr"
        },
        "debug_file": { # doesn't include debug messages from imported modules
                 "class": "logging.FileHandler",
                 "level": "DBG_SCR",
                 "formatter": "file",
                 "filename": f"{args.output_prefix}.debug_file.log",
                 "mode": "w"
                },
        "low_lvl_debug_file": { # doesn't include debug messages from imported modules
                 "class": "logging.FileHandler",
                 "level": "LOW_DBG_SCR",
                 "formatter": "file",
                 "filename": f"{args.output_prefix}.low_lvl_debug_file.log",
                 "mode": "w"
                },
        "full_file": { #includes debug messages from imported modules
            "class": "logging.FileHandler",
            "formatter": "file",
            "filename": f"{args.output_prefix}.total.log",
            "mode": "w"
        }
    },
    "root": {
        "level": "DEBUG",
        "handlers": [
            "stderr",
            "stdout",
            "debug_file",
            "low_lvl_debug_file",
            "full_file"
        ]
    }
}

def filter_maker(level):
    level = getattr(logging, level)

    def filter(record):
        return record.levelno <= level

    return filter

logging.config.dictConfig(log_config_dict)
# -------- End of Config loggers --------

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

syn_file_key_column, syn_file_value_column = args.syn_file_key_column, args.syn_file_value_column

#label_dict = {genome: label for genome, label in zip(args.genome_orderlist,
#                                                     args.genome_orderlist if args.genome_labellist is None else args.genome_labellist)}

synteny_format = args.synteny_format

# read files
logging.info(f"Data directory: {data_dir_path}/\n")
#print(data_dir_path)
# whitelist is obligatory

logging.info(f"Genomes:")
for genome in genome_orderlist:
    logging.info(TAB * TAB_NUM_PER_LVL + f'{genome}')
logging.info("\n")

figure_config = {}

default_config = {"scaffold": {
                               "show":               True,
                               "rounded":            True,
                               "stranded":           False,
                               "stranded_end":       False,
                               "show_label":         True,
                               "label_font":         None,
                               "label_fontsize":     7 if args.scaffold_label_fontsize is None else args.scaffold_label_fontsize,
                               "label_fontweight":   None,
                               "label_angle":        0 if args.scaffold_label_angle is None else args.scaffold_label_angle,
                               "fill_color":         None,
                               "fill":               True,
                               "edge":               False,
                               "edge_color":         None,
                               "middle_line":        True,
                               "middle_line_color":  None,
                               },
                  "superscaffold": {
                                    "show": False, # Superscaffold is not shown by default
                                    "rounded": True,
                                    "stranded": False,
                                    "stranded_end": False,
                                    "show_label": True,
                                    "label_font": None,
                                    "label_fontsize": 7,
                                    "label_fontweight": None,
                                    "label_angle": 0,
                                    "fill_color": None,
                                    "fill": True,
                                    "edge": False,
                                    "edge_color": None,
                                    "middle_line": True,
                                    "middle_line_color": None,
                                    },
                  "genome": {
                             "show_label":         True,
                             "label_font":         None,
                             "label_fontsize":     11 if args.genome_label_fontsize is None else args.genome_label_fontsize,
                             "label_fontweight":   None,
                             "label_angle":        0 if args.genome_label_angle is None else args.genome_label_angle,
                             "fill_color":         None, #"#10d541" # green
                             "edge_color":         None
                            },
                  "coordinate_system": "scaffold"
                  }

number_of_genomes = len(genome_orderlist)
color_list = distinctipy.get_colors(number_of_genomes)

default_genome_color_dict = {genome_orderlist[genome_index]: rgb_tuple_to_hex(color_list[genome_index]) for genome_index in range(0, number_of_genomes)}

#if genome_colors:
#    for g_color_index in range(0, len(genome_colors)):
#        if genome_colors[g_color_index] is not None:
#            colors[g_color_index] = genome_colors[g_color_index]

#color_list = list(map(rgb_tuple_to_hex, colors))
#color_df = pd.DataFrame.from_dict({"genome_id": genome_orderlist,
#                                   "color": color_list},)
#color_df.set_index("genome_id", inplace=True)
#color_df.to_csv(f"{args.output_prefix}.genome.colors", sep="\t", index=True, header=True)

config_dict = {genome: {} for genome in genome_orderlist}

for genome_index in range(0, number_of_genomes):
    genome = genome_orderlist[genome_index]
    config_dict[genome] = {
                           "scaffold":      {"general": {},
                                             "syn_df": None,
                                             "settings_df": None,
                                             "len_df": None,
                                             "centromere_df": None,
                                             "color_df": None,
                                             "whitelist": [],
                                             "invertlist": [],
                                             "orderlist": [],
                                             "queryswitchstrandlist": [],
                                             "targetswitchstrandlist": []
                                             },
                           "superscaffold": {"general": {},
                                             "syn_df": None,
                                             "settings_df": None,
                                             "len_df": None,
                                             "centromere_df": None,
                                             "color_df": None,
                                             "whitelist": [],
                                             "invertlist": [],
                                             "orderlist": [],
                                             "queryswitchstrandlist": [],
                                             "targetswitchstrandlist": []
                                             },
                           "genome": {},
                           "coordinate_system": None
                          }

extension_dict = {
                  "general_config_file": ["general.config.yaml"],
                  "superscaffold_file": ["agp"],
                  "config_file": ["config.tsv"],
                  "whitelist": ["whitelist"],
                  "orderlist": ["orderlist"],
                  "invertlist": ["invertlist"],
                  "queryswitchstrandlist": ["queryswitchstrandlist"],
                  "targetswitchstrandlist": ["targetswitchstrandlist"],
                  "syn": ["syn"],
                  "len": ["len"],
                  "colors": ["colors"],
                  "centromere": ["centromere.bed"]
                  }


def add_datatype_to_ext(datatype_f, extension_list):
    return list(map(lambda ext_f: f"{datatype_f}.{ext_f}", extension_list))


superscaffold_df_dict = {}
logging.info(f"Parsing input files...")
logging.dbg_scr(f"Input files:")

for genome_index in range(0, number_of_genomes):
    genome = genome_orderlist[genome_index]
    # nonempty whitelist file is necessary for each genome
    logging.dbg_scr(TAB * TAB_NUM_PER_LVL + f'{genome}:')

    # general config file might be absent, but should be parsable by yaml module
    general_config_file = get_filenames_for_extension(data_dir_path / genome, extension_list=extension_dict["general_config_file"])
    if general_config_file is None:
        logging.dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + f'general config: ABSENT (no )')
    else:
        logging.dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + f'general config: {general_config_file}')
        with open(general_config_file, 'r') as general_config_fd:
            tmp_dict = yaml.safe_load(general_config_fd)
        if "scaffold" in tmp_dict:
            config_dict[genome]["scaffold"]["general"] = tmp_dict["scaffold"]
        if "superscaffold" in tmp_dict:
            config_dict[genome]["superscaffold"]["general"] = tmp_dict["superscaffold"]
        if "genome" in tmp_dict:
            config_dict[genome]["genome"] = tmp_dict["genome"]
        if "coordinate_system" in tmp_dict:
            config_dict[genome]["coordinate_system"] = tmp_dict["coordinate_system"]

    # superscaffold file (agp) might be absent or empty
    # TODO: add parsing for other formats
    superscaffold_file = get_filenames_for_extension(data_dir_path / genome, extension_list=extension_dict["superscaffold_file"])
    if superscaffold_file is None:
        logging.dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + f'superscaffold file: ABSENT')
        superscaffold_df_dict[genome] = None
    else:
        try:
            superscaffold_df_dict[genome] = CollectionAGP(in_file=superscaffold_file)
            logging.dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + f'superscaffold file: {superscaffold_file}')
            # logging.warn_scr(f"AGP file is present for {genome}!!! Ignoring syn, whitelist, invertlist, orderlist, targetswitchstrandlist and queryswitchstrandlist for this genome! "
            #                 f"Seeking for agp_syn, agp_whitelist, agp_invertlist, agp_orderlist, agp_targetswitchstrandlist, agp_queryswitchstrandlist and agp_colors. "
            #                 f"IMPORTANT!!! Len list for contigs is still necessary.\n")
        except pd.errors.EmptyDataError:
            logging.dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + f'superscaffold file: EMPTY')
            # raise pd.errors.EmptyDataError(f"ERROR!!! AGP file for {genome} is present, but empty! Replace it by the valid agp file or remove it!")

    for datatype in "superscaffold", "scaffold":
        logging.dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + f'datatype: {datatype}')

        # datatype config file might be absent or empty or parsable into pandas dataframe with first row being a header
        datatype_config_file = get_filenames_for_extension(data_dir_path / genome, extension_list=add_datatype_to_ext(datatype, extension_dict["config_file"]))
        if datatype_config_file is None:
            logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'config: ABSENT')
        else:
            try:
                config_dict[genome][datatype]["settings_df"] = pd.read_csv(datatype_config_file, sep="\t", header=0, comment=None, index_col=0)
                logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'config: {datatype_config_file}')
            except pd.errors.EmptyDataError:
                config_dict[genome][datatype]["settings_df"] = None
                logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'config: EMPTY')

        #  list files below might be abent or empty
        for list_type in "whitelist", "orderlist", "invertlist", "queryswitchstrandlist", "targetswitchstrandlist":
            input_file = get_filenames_for_extension(data_dir_path / genome, extension_list=add_datatype_to_ext(datatype, extension_dict[list_type]))
            if input_file is None:
                logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'{list_type}: ABSENT')
                config_dict[genome][datatype][list_type] = pd.Series(dtype=str)
            else:
                try:
                    config_dict[genome][datatype][list_type] = pd.read_csv(input_file, sep="\t", header=None, comment="#").squeeze("columns")
                    logging.dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + f'{list_type}: {input_file}')
                except pd.errors.EmptyDataError:
                    logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'{list_type}: EMPTY')
                    config_dict[genome][datatype][list_type] = pd.Series(dtype=str)

        # syn file might be empty or absent
        syn_file = get_filenames_for_extension(data_dir_path / genome, extension_list=add_datatype_to_ext(datatype, extension_dict["syn"]))
        if syn_file is None:
            logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'syn: ABSENT')
            config_dict[genome][datatype]["syn_df"] = pd.DataFrame(columns=["syn"])
        else:
            try:
                config_dict[genome][datatype]["syn_df"] = pd.read_csv(syn_file, usecols=(syn_file_key_column, syn_file_value_column),
                                                                      sep="\t", header=None, comment="#",
                                                                      names=["key", "syn"] if syn_file_key_column <= syn_file_value_column else ["syn", "key"]).set_index("key")
                logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'syn: {syn_file}')
            except pd.errors.EmptyDataError:
                logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'syn: EMPTY')
                config_dict[genome][datatype]["syn_df"] = pd.DataFrame(columns=["syn"])
        config_dict[genome][datatype]["renamed_whitelist"] = config_dict[genome][datatype]["whitelist"].replace(config_dict[genome][datatype]["syn_df"]["syn"].to_dict(), inplace=False)
        # centromere.bed might be empty or absent
        centromere_file = get_filenames_for_extension(data_dir_path / genome, extension_list=add_datatype_to_ext(datatype, extension_dict["centromere"]))
        if centromere_file is None:
            config_dict[genome][datatype]["centromere_df"] = pd.DataFrame(columns=["scaffold_id", "start", "end"])
            config_dict[genome][datatype]["centromere_df"].set_index("scaffold_id", inplace=True)
            logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'centromeres: ABSENT')
        else:
            try:
                config_dict[genome][datatype]["centromere_df"] = pd.read_csv(centromere_file, usecols=(0, 1, 2), index_col=0,
                                                                             header=None, sep="\t", names=["scaffold_id", "start", "end"])
                logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'centromeres: {centromere_file}')
            except pd.errors.EmptyDataError:
                logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'centromeres: EMPTY')
                config_dict[genome][datatype]["centromere_df"] = pd.DataFrame(columns=["scaffold_id", "start", "end"])
                config_dict[genome][datatype]["centromere_df"].set_index("scaffold_id", inplace=True)
        config_dict[genome][datatype]["centromere_df"].rename(index=config_dict[genome][datatype]["syn_df"]["syn"].to_dict(), inplace=True)

        # color file might be empty or absent
        color_file = get_filenames_for_extension(data_dir_path / genome, extension_list=add_datatype_to_ext(datatype, extension_dict["colors"]))
        if color_file is None:
            logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'colors: ABSENT')
            config_dict[genome][datatype]["color_df"] = pd.DataFrame(columns=["color"])
        else:
            try:
                config_dict[genome][datatype]["color_df"] = pd.read_csv(color_file,sep="\t", header=0, names=["scaffold", "color"], index_col=0)
                logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'colors: {color_file}')
            except pd.errors.EmptyDataError:
                logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'colors: EMPTY')
                config_dict[genome][datatype]["color_df"] = pd.DataFrame(columns=["scaffold_id", "color"])
                config_dict[genome][datatype]["color_df"].set_index("scaffold_id", inplace=True)

    # nonempty lenlist file is necessary for each genome
    lenlist_file = get_filenames_for_extension(data_dir_path / genome, extension_list=add_datatype_to_ext("scaffold", extension_dict["len"]))
    try:
        config_dict[genome]["scaffold"]["len_df"] = pd.read_csv(lenlist_file, sep="\t", header=None, comment="#",
                                                                       names=["scaffold", "length"], index_col=0)
        logging.dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + f'len: {lenlist_file}')
    except pd.errors.EmptyDataError:
        logging.dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f'len: EMPTY')
        raise pd.errors.EmptyDataError(f"ERROR!!! lenlist for {genome} is empty. add scaffold ids and its lengths to it!")
    except ValueError:
        raise ValueError(f"Lenlist for {genome} is absent!")

    # filter and rename scaffolds in superscaffold_file
    if superscaffold_df_dict[genome] is not None:
        superscaffold_df_dict[genome].filter_records(white_list=config_dict[genome][datatype]["whitelist"])
        config_dict[genome]["superscaffold"]["len_df"] = superscaffold_df_dict[genome].get_seq_len()
        superscaffold_df_dict[genome] = superscaffold_df_dict[genome].get_nongap_records().reset_index(drop=True).set_index("scaffold")
        superscaffold_df_dict[genome].rename(index=config_dict[genome][datatype]["syn_df"]["syn"].to_dict(), inplace=True)

        config_dict[genome]["scaffold"]["whitelist"] = superscaffold_df_dict[genome]["part_id/gap_length"]
        config_dict[genome]["scaffold"]["orderlist"] = pd.concat([superscaffold_df_dict[genome].loc[agp_scaffold_id,
                                                                                                    "part_id/gap_length"] for agp_scaffold_id in config_dict[genome]["superscaffold"]["orderlist"]])
        config_dict[genome]["scaffold"]["invertlist"] = superscaffold_df_dict[genome].loc[superscaffold_df_dict[genome]["orientation/evidence"] == "-", "part_id/gap_length"]
        config_dict[genome]["scaffold"]["queryswitchstrandlist"] = superscaffold_df_dict[genome].loc[config_dict[genome]["superscaffold"]["queryswitchstrandlist"],
                                                                                                    "part_id/gap_length"]
        config_dict[genome]["scaffold"]["targetswitchstrandlist"] = superscaffold_df_dict[genome].loc[config_dict[genome]["superscaffold"]["targetswitchstrandlist"],
                                                                                                     "part_id/gap_length"]
    else:
        for metadata_type in ["len_df", "whitelist", "orderlist", "invertlist", "queryswitchstrandlist", "queryswitchstrandlist"]:
            config_dict[genome]["superscaffold"][metadata_type] = None

logging.dbg_scr("\n")
# -------

logging.info("Resolving settings...")

if args.settings_level == 'default':
    logging.info(TAB * 1 * TAB_NUM_PER_LVL + "Settings level was set to default! Ignoring genome and scaffold settings")
    for genome in genome_orderlist:
        for datatype in "superscaffold", "scaffold":
            config_dict[genome][datatype]["general"] = default_config[datatype]

            config_dict[genome][datatype]["settings_df"] = pd.DataFrame(index=config_dict[genome][datatype]["renamed_whitelist"],)
            for parameter in config_dict[genome][datatype]["general"]:
                config_dict[genome][datatype]["settings_df"][parameter] = config_dict[genome][datatype]["general"][parameter]

        config_dict[genome]["genome"] = default_config["genome"]
        config_dict[genome]["genome"]["fill_color"] = default_genome_color_dict[genome]
        config_dict[genome]["genome"]["edge_color"] = default_genome_color_dict[genome]
        config_dict[genome]["coordinate_system"] = default_config["coordinate_system"]

else:
    # filling absent settings in general configs from default settings
    for genome in genome_orderlist:
        for datatype in "superscaffold", "scaffold":
            copy_absent_entries(default_config[datatype], config_dict[genome][datatype]["general"])
        if config_dict[genome]['coordinate_system'] is None:
            config_dict[genome]["coordinate_system"] = default_config["coordinate_system"]
        if ("fill_color" not in config_dict[genome]["genome"]) or (config_dict[genome]["genome"]["fill_color"] is None):
            config_dict[genome]["genome"]["fill_color"] = default_genome_color_dict[genome]
        if ("edge_color" not in config_dict[genome]["genome"]) or (config_dict[genome]["genome"]["fill_color"] is None):
            config_dict[genome]["genome"]["edge_color"] = default_genome_color_dict[genome]
        copy_absent_entries(default_config["genome"], config_dict[genome]["genome"])

    if args.settings_level == 'genome':
        for genome in genome_orderlist:
            for datatype in "superscaffold", "scaffold":
                config_dict[genome][datatype]["settings_df"] = pd.DataFrame(index=config_dict[genome][datatype]["renamed_whitelist"],)
                for parameter in config_dict[genome][datatype]["general"]:
                    config_dict[genome][datatype]["settings_df"][parameter] = config_dict[genome][datatype]["general"][parameter]

                config_dict[genome][datatype]["settings_df"]["fill_color"] = config_dict[genome]["genome"]["fill_color"]
                config_dict[genome][datatype]["settings_df"]["edge_color"] = config_dict[genome]["genome"]["edge_color"]
    elif args.settings_level == 'scaffold':
        for genome in genome_orderlist:
            #print("AAAa")
            for datatype in "superscaffold", "scaffold":
                if config_dict[genome][datatype]["settings_df"] is None:
                    config_dict[genome][datatype]["settings_df"] = pd.DataFrame(index=config_dict[genome][datatype]["renamed_whitelist"],)
                    for parameter in config_dict[genome][datatype]["general"]:
                        config_dict[genome][datatype]["settings_df"][parameter] = config_dict[genome][datatype]["general"][parameter]

                else:
                    #print(config_dict[genome][datatype]["settings_df"])
                    whitelist_set = set(config_dict[genome][datatype]["renamed_whitelist"])
                    initial_entries_set = set(config_dict[genome][datatype]["settings_df"].index)
                    present_entries_set = initial_entries_set & whitelist_set
                    absent_entries_list = list(whitelist_set - present_entries_set)

                    common_parameters_set = set(config_dict[genome][datatype]["general"].keys())
                    present_parameters_set = set(config_dict[genome][datatype]["settings_df"].columns)
                    absent_parameters_set = common_parameters_set - present_parameters_set

                    # keep only whitelisted scaffolds/superscaffolds and registered parameters.
                    config_dict[genome][datatype]["settings_df"] = config_dict[genome][datatype]["settings_df"].loc[config_dict[genome][datatype]["settings_df"].index.isin(config_dict[genome][datatype]["renamed_whitelist"]),
                                                                                                                    config_dict[genome][datatype]["settings_df"].columns.isin(config_dict[genome][datatype]["general"].keys())]
                    for parameter in absent_parameters_set:  # add missing parameters to the existing entries
                        #print()
                        config_dict[genome][datatype]["settings_df"].loc[:, parameter] = config_dict[genome][datatype]["general"][parameter]
                    if absent_entries_list:
                        #print(config_dict[genome][datatype]["settings_df"].columns)
                        absent_entries_empty_df = pd.DataFrame(index=absent_entries_list, columns=config_dict[genome][datatype]["settings_df"].columns)
                        #print(absent_entries_empty_df)
                        config_dict[genome][datatype]["settings_df"] = pd.concat((config_dict[genome][datatype]["settings_df"], absent_entries_empty_df), axis=0)
                        #print(config_dict[genome][datatype]["settings_df"])
                        for parameter in config_dict[genome][datatype]["general"].keys():  # add missing entries
                            config_dict[genome][datatype]["settings_df"].loc[absent_entries_list, parameter] = config_dict[genome][datatype]["general"][parameter]
                #print(config_dict[genome][datatype]["settings_df"]["fill_color"].isnull())
                print(config_dict[genome][datatype]["settings_df"])
                #print(list(config_dict[genome][datatype]["settings_df"].columns))

                #print(config_dict[genome]["genome"]["fill_color"])
                #print(config_dict[genome]["genome"])
                #print(config_dict[genome][datatype]["settings_df"]["fill_color"].isnull())
                config_dict[genome][datatype]["settings_df"].loc[config_dict[genome][datatype]["settings_df"]["fill_color"].isnull(), "fill_color"] = config_dict[genome]["genome"]["fill_color"]
                print(config_dict[genome][datatype]["settings_df"])
                # set edge_color same to fill_color if it is not set
                config_dict[genome][datatype]["settings_df"].loc[config_dict[genome][datatype]["settings_df"]["edge_color"].isnull(), "edge_color"] = config_dict[genome][datatype]["settings_df"].loc[config_dict[genome][datatype]["settings_df"]["edge_color"].isnull(), "fill_color"]
                print(config_dict[genome][datatype]["settings_df"])
                #if config_dict[genome][datatype]["color_df"] =
    else:
        raise ValueError(f"ERROR!!! Unknown settings level (f{args.settings_level})")

logging.low_dbg_scr("Printing general config files...")

for genome in genome_orderlist:
    logging.low_dbg_scr(TAB * 1 * TAB_NUM_PER_LVL + f"{genome}:")
    for datatype in "scaffold", "superscaffold", "genome":
        logging.low_dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + f"{datatype}:\n")
        logging.low_dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + "\n" + str(config_dict[genome][datatype]) + "\n")
    logging.low_dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + f"coordinate_system: {config_dict[genome]['coordinate_system']}\n")

for genome in superscaffold_df_dict:
    logging.low_dbg_scr(TAB * 1 * TAB_NUM_PER_LVL + f"{genome}:")
    if superscaffold_df_dict[genome] is not None:
        logging.low_dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + "\n" + str(superscaffold_df_dict[genome]))
    else:
        logging.low_dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + "ABSENT")

logging.low_dbg_scr("Printing agp files...")
for genome in superscaffold_df_dict:
    logging.low_dbg_scr(TAB * 1 * TAB_NUM_PER_LVL + f"{genome}:")
    if superscaffold_df_dict[genome] is not None:
        logging.low_dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + "\n" + str(superscaffold_df_dict[genome]))
    else:
        logging.low_dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + "ABSENT")
        # ---------------------------- Read genome config if it was set --------------------------------------
logging.low_dbg_scr("\n")

logging.low_dbg_scr("Printing config dict...")
for genome in genome_orderlist:
    logging.low_dbg_scr(TAB * 1 * TAB_NUM_PER_LVL + f"{genome}:")
    for metadata_type in ["len_df", "whitelist", "orderlist", "invertlist", "queryswitchstrandlist", "targetswitchstrandlist"]:
        logging.low_dbg_scr(TAB * 2 * TAB_NUM_PER_LVL + f"{metadata_type}:")
        for datatype in ["scaffold", "superscaffold"]:
            logging.low_dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f"{datatype}:")
            #print(config_dict[genome][datatype])
            if config_dict[genome][datatype][metadata_type] is None:
                logging.low_dbg_scr(TAB * 4 * TAB_NUM_PER_LVL + "ABSENT")
            else:
                logging.low_dbg_scr(TAB * 4 * TAB_NUM_PER_LVL + "\n" + str(config_dict[genome][datatype][metadata_type]))
        logging.low_dbg_scr(TAB * 3 * TAB_NUM_PER_LVL + f"coordinate system: {config_dict[genome]['coordinate_system']}")


def check_for_inversion_and_strand_switch_request(scaffold_id, strand_switch_symbol, inversion_symbol,):
    request_dict = {"inversion": False,
                    "query_switch_strand": False,
                    "target_switch_strand": False}
    processed_scaffold_id = scaffold_id

    if scaffold_id[-1] == inversion_symbol:
        request_dict["inversion"] = True
        processed_scaffold_id = scaffold_id[:-1]
    if processed_scaffold_id[-3:] == strand_switch_symbol * 3:
        request_dict["query_switch_strand"] = True
        request_dict["target_switch_strand"] = True
        processed_scaffold_id = processed_scaffold_id[:-3]
    elif processed_scaffold_id[-2:] == strand_switch_symbol * 2:
        request_dict["target_switch_strand"] = True
        processed_scaffold_id = processed_scaffold_id[:-2]
    elif processed_scaffold_id[-1:] == strand_switch_symbol * 1:
        request_dict["query_switch_strand"] = True
        processed_scaffold_id = processed_scaffold_id[:-1]

    if len(processed_scaffold_id) == 0:
        raise ValueError(f"ERROR!!! No symbols were left after processing of the scaffold id {scaffold_id} !")

    return processed_scaffold_id, request_dict


def get_original_scaffold_id(syn_df, scaffold_id):
    if scaffold_id not in list(syn_df["syn"]):
        #print(syn_df["syn"])
        return scaffold_id
    return syn_df.index[list(syn_df["syn"]).index(scaffold_id)]

# count number of chromosomes/scaffold and total length of them
# Check if adjustment for coordinate system is necessary
logging.info("Coordinate systems:")
for genome in genome_orderlist:
    logging.low_dbg_scr(TAB * 1 * TAB_NUM_PER_LVL + f"{genome}:")
    if ("show" in config_dict[genome]["superscaffold"]["general"]) and config_dict[genome]["superscaffold"]["general"]["show"]:
        if config_dict[genome]["coordinate_system"] != "superscaffold":
            config_dict[genome]["coordinate_system"] = "superscaffold"
            logging.warn_scr(TAB * 2 * TAB_NUM_PER_LVL + f"ENFORCING superscaffold coordinate system as superscaffolds were requested to be shown for this genome")
        else:
            logging.info(TAB * 2 * TAB_NUM_PER_LVL + f"{config_dict[genome]['coordinate_system']}")
    else:
        logging.info(TAB * 2 * TAB_NUM_PER_LVL + f"{config_dict[genome]['coordinate_system']}")

if args.genome_config: # TODO: AFTRER RECENT CHANGES GENOME CONFIG NO MORE WORKS CORRECTLY. FIX IT
    logging.warn_scr("Genome config file is set!!! Ignoring invertlist, orderlist, targetswitchstrandlist and queryswitchstrandlist for all genomes!\n")

    genome_config_tmp_dict = OrderedDict()

    with open(args.genome_config, "r") as in_fd:
        for line in in_fd:
            if line == "\n":
                continue
            line_list = line.strip().split("\t", 2)
            genome_colors.append(line_list[1] if line_list[1] != "." else None) # read_colors
            genome_config_tmp_dict[line_list[0]] = list(map(lambda s: s.split(","), line_list[2].split("\t")))

    for genome in genome_config_tmp_dict:
        datatype = config_dict[genome]["coordinate_system"]

        config_dict[genome][datatype]["invertlist"] = []
        config_dict[genome][datatype]["orderlist"] = []
        config_dict[genome][datatype]["targetswitchstrandlist"] = []
        config_dict[genome][datatype]["queryswitchstrandlist"] = []
        for entry in genome_config_tmp_dict[genome]:
            #print(entry)
            for scaffold_id in entry:
                #print(scaffold_id)

                processed_scaffold_id, request_dict = check_for_inversion_and_strand_switch_request(scaffold_id,
                                                                                                    args.strand_switch_label,
                                                                                                    args.inverted_scaffold_label)
                #print(processed_scaffold_id, request_dict)
                #print(get_original_scaffold_id(syn_df_dict[genome], processed_scaffold_id))
                config_dict[genome][datatype]["orderlist"].append(processed_scaffold_id)
                if request_dict["inversion"]:
                    config_dict[genome][datatype]["invertlist"].append(processed_scaffold_id)

                if request_dict["query_switch_strand"]:
                    config_dict[genome][datatype]["queryswitchstrandlist"].append(get_original_scaffold_id(config_dict[genome][datatype]["syn_df"],
                                                                                               processed_scaffold_id))
                if request_dict["target_switch_strand"]:
                    config_dict[genome][datatype]["targetswitchstrandlist"].append(get_original_scaffold_id(config_dict[genome][datatype]["syn_df"],
                                                                                                processed_scaffold_id))

        config_dict[genome][datatype]["invertlist"] = pd.Series(config_dict[genome][datatype]["invertlist"], dtype='str')
        config_dict[genome][datatype]["orderlist"] = pd.Series(config_dict[genome][datatype]["orderlist"], dtype='str')
        config_dict[genome][datatype]["queryswitchstrandlist"] = pd.Series(config_dict[genome][datatype]["queryswitchstrandlist"], dtype='str')
        config_dict[genome][datatype]["targetswitchstrandlist"] = pd.Series(config_dict[genome][datatype]["targetswitchstrandlist"], dtype='str')

#filter len list
for genome in genome_orderlist:
    config_dict[genome]["scaffold"]["len_df"] = config_dict[genome]["scaffold"]["len_df"].loc[config_dict[genome]["scaffold"]["len_df"].index.isin(config_dict[genome]["scaffold"]["whitelist"])]

# rename len lists
for genome in genome_orderlist:
    for datatype in "scaffold", "superscaffold":
        if not config_dict[genome][datatype]["syn_df"].empty:
            config_dict[genome][datatype]["len_df"].rename(index=config_dict[genome][datatype]["syn_df"]["syn"].to_dict(), inplace=True)

for genome in genome_orderlist:
    for datatype in "scaffold", "superscaffold":
        if config_dict[genome][datatype]["len_df"] is not None:
            config_dict[genome][datatype]["len_df"] = config_dict[genome][datatype]["len_df"].reindex(config_dict[genome][datatype]["orderlist"]).dropna()

total_len_dict = {genome: sum(config_dict[genome][config_dict[genome]["coordinate_system"]]["len_df"]["length"]) for genome in genome_orderlist}
chr_number_dict = {genome: len(config_dict[genome][config_dict[genome]["coordinate_system"]]["len_df"]) for genome in genome_orderlist}

max_genome_length = max(list(total_len_dict.values()))

if synteny_format == "psl":
    logging.info("Synteny file format: psl\n")
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
        logging.info(TAB * 1 * TAB_NUM_PER_LVL + f"Query genome: {genome_orderlist[genome_index]}")
        logging.info(TAB * 1 * TAB_NUM_PER_LVL + "Target genome: {0}".format(genome_orderlist[genome_index + 1]))
        logging.info(TAB * 1 * TAB_NUM_PER_LVL + "psl file: {0}".format(get_filenames_for_extension(data_dir_path / genome_orderlist[genome_index],
                                                                                                    extension_list=["psl", "psl.gz"])))
        logging.info(TAB * 1 * TAB_NUM_PER_LVL + "Parsing...")

        synteny_dict[genome_orderlist[genome_index]] = CollectionPSL(get_filenames_for_extension(data_dir_path / genome_orderlist[genome_index],
                                                                                                 extension_list=["psl", "psl.gz"]),
                                                                     target_white_list=config_dict[genome_orderlist[genome_index + 1]]["scaffold"]["whitelist"],
                                                                     query_white_list=config_dict[genome_orderlist[genome_index]]["scaffold"]["whitelist"],
                                                                     #target_syn_dict=syn_df_dict[genome_orderlist[genome_index + 1]]["syn"].to_dict(),
                                                                     #query_syn_dict=syn_df_dict[genome_orderlist[genome_index]]["syn"].to_dict(),
                                                                     parsing_mode="coordinates_only").records.sort_values(by=[query_scaffold_id_column_name,
                                                                                                                              query_start_column_name,
                                                                                                                              query_end_column_name,
                                                                                                                              target_scaffold_id_column_name,
                                                                                                                              target_start_column_name,
                                                                                                                              target_end_column_name])
        logging.info("\n")
else:
    raise ValueError("ERROR!!! {0} format is not implemented yet!".format("psl"))

#----------------------------- Assignment of IDs to synteny blocks -------------------------------
for genome in synteny_dict:
    synteny_dict[genome]['synteny_block_id'] = ["SB_{0}".format(block_id) for block_id in range(1, len(synteny_dict[genome]) + 1)]
    #print(synteny_dict[genome])
#-------------------------------------------------------------------------------------------------

#print(syn_df_dict)
#print(synteny_dict["dingo"])
genome_number = len(genome_orderlist)

#------------------ Preprocessing of synteny blocks -----------------------------------
##----------------------- Detect nested blocks ----------------------------------------
logging.info("Preprocessing of synteny blocks\n")


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
    #print(tmp_dict)
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
    #print(tmp_dict)
    #print("AAAAAAAAAA")
    #print(synteny_dict[genome])
    if args.remove_same_coords_blocks and not tmp_dict[genome].empty:
        block_remove_dict[genome] = list(tmp_dict[genome]["target_same_coords"].dropna())

    if args.remove_nested_blocks and not tmp_dict[genome].empty:
        block_remove_dict[genome] += list(tmp_dict[genome][tmp_dict[genome]["query_nested_in"].notna()]["synteny_block_id"])
        block_remove_dict[genome] += list(tmp_dict[genome][tmp_dict[genome]["target_nested_in"].notna()]["synteny_block_id"])

    block_remove_dict[genome] = set(block_remove_dict[genome])
    synteny_dict[genome] = synteny_dict[genome][~synteny_dict[genome]["synteny_block_id"].isin(block_remove_dict[genome])]


#print(block_remove_dict)
#------------------ Classification of synteny blocks ----------------------------------
logging.info("Classification of synteny blocks...\n")
for genome_index in range(0, genome_number - 1): # genome_orderlist[:-1]:  # all query genomes
    target_genome_index = genome_index + 1
    genome = genome_orderlist[genome_index]
    target_genome = genome_orderlist[target_genome_index]
    logging.info(TAB * 1 * TAB_NUM_PER_LVL + f"Query: {genome}")
    logging.info(TAB * 1 * TAB_NUM_PER_LVL + f"Target: {target_genome}")
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

        logging.info(TAB * 2 * TAB_NUM_PER_LVL + "Flag --invert_major_strand flag was set. Switching major strand for all genomes and all chromosomes...")
        major_strand_series[minus_strand_index] = "+"
        major_strand_series[plus_strand_index] = "-"

    # apply a local major strand switch

    if not config_dict[genome]["scaffold"]["queryswitchstrandlist"].empty:
        logging.info(TAB * 2 * TAB_NUM_PER_LVL + f"Switching major strand for {genome} as query for {genome} vs {target_genome} alignment...")
        #print(major_strand_series)
        #switch_strand_scaffolds = set(switchstrandlist_series_dict[genome]) & set(major_strand_series.index.get_level_values(0))
        #print(major_strand_series == "-")
        switch_series = pd.Series(major_strand_series.index.get_level_values(query_scaffold_id_column_name).isin(config_dict[genome]["scaffold"]["queryswitchstrandlist"]))
        switch_series.index = major_strand_series.index

        minus_strand_index = switch_series & (major_strand_series == "-")
        plus_strand_index = switch_series & (major_strand_series == "+")
        major_strand_series[minus_strand_index] = "+"
        major_strand_series[plus_strand_index] = "-"

    # apply switchstrand list of target genome.
    #print(metadata_dict["scaffold"]["targetswitchstrandlist_series_dict"][target_genome])
    if not config_dict[target_genome]["scaffold"]["targetswitchstrandlist"].empty:
        logging.info(TAB * 2 * TAB_NUM_PER_LVL + f"Switching major strand for {target_genome} as target for {genome} vs {target_genome} alignment...")
        switch_series = pd.Series(major_strand_series.index.get_level_values(target_scaffold_id_column_name).isin(config_dict[target_genome]["scaffold"]["targetswitchstrandlist"]))
        switch_series.index = major_strand_series.index

        minus_strand_index = switch_series & (major_strand_series == "-")
        plus_strand_index = switch_series & (major_strand_series == "+")
        major_strand_series[minus_strand_index] = "+"
        major_strand_series[plus_strand_index] = "-"

    synteny_dict[genome]["major_strand"] = major_strand_series
    #print("CCCC")
    #print(major_strand_series)
    synteny_dict[genome].reset_index(level=1, drop=False, inplace=True)
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
    genome = genome_orderlist[genome_index]
    if superscaffold_df_dict[genome] is None:
        synteny_dict[genome_orderlist[genome_index]][target_scaffold_id_column_name] = synteny_dict[genome_orderlist[genome_index]][target_scaffold_id_column_name].replace(config_dict[genome_orderlist[genome_index + 1]]["scaffold"]["syn_df"]["syn"].to_dict())
        synteny_dict[genome_orderlist[genome_index]][query_scaffold_id_column_name] = synteny_dict[genome_orderlist[genome_index]][query_scaffold_id_column_name].replace(config_dict[genome_orderlist[genome_index]]["scaffold"]["syn_df"]["syn"].to_dict())
    synteny_dict[genome_orderlist[genome_index]] = synteny_dict[genome_orderlist[genome_index]].sort_values(by=[query_scaffold_id_column_name,
                                                                                                                query_start_column_name,
                                                                                                                query_end_column_name,
                                                                                                                target_scaffold_id_column_name,
                                                                                                                target_start_column_name,
                                                                                                                target_end_column_name])
    for datatype in "scaffold", "superscaffold":
        config_dict[genome][datatype]["settings_df"].rename(index=config_dict[genome][datatype]["syn_df"]["syn"].to_dict(),
                                                            inplace=True)
    #for index, genome in zip(range(0, len(genome_orderlist) - 1), genome_orderlist[:-1]):
#    synteny_dict[genome].to_csv("{0}.{1}.to.{2}.raw.renamed.tab".format(args.output_prefix,
#                                                                        genome,
#                                                                        genome_orderlist[index + 1]),
#                                sep="\t", index=False, header=True)

#--------------------------------------------------------------------------------------

if args.remove_scaffolds_absent_in_orderlist:
    for index in range(0, len(genome_orderlist) - 1):
        query_genome = genome_orderlist[index]
        target_genome = genome_orderlist[index + 1]
        bool_series = synteny_dict[query_genome][query_scaffold_id_column_name].isin(config_dict[query_genome]["scaffold"]["orderlist"])
        #print(sum(bool_series))
        bool_series &= synteny_dict[query_genome][target_scaffold_id_column_name].isin(config_dict[target_genome]["scaffold"]["orderlist"])
        #print(sum(bool_series))
        synteny_dict[query_genome] = synteny_dict[query_genome][bool_series]

#-------------------------- Inversion of coordinates ----------------------------------
logging.info("Inverting scaffolds (if necessary)...")
for genome_index in range(0, genome_number):
    genome = genome_orderlist[genome_index]
    logging.info(TAB * 1 * TAB_NUM_PER_LVL + genome + ":")

    #print("Inverting (if necessary) {0} scaffolds...".format(genome))
    #print(centromere_df_dict[genome])
    #print(config_dict[genome]["scaffold"]["len_df"])

    if genome_index < (genome_number - 1):  # apply listed inversion for all query genomes, i.e. for all except the last genome
        logging.info(TAB * 2 * TAB_NUM_PER_LVL + "Inverting query coordinates in synteny file...")
        #print(config_dict[genome]["scaffold"]["len_df"])
        synteny_dict[genome] = invert_coordinates_in_synteny_table(synteny_dict[genome],
                                                                   config_dict[genome]["scaffold"]["invertlist"],
                                                                   config_dict[genome]["scaffold"]["len_df"],
                                                                   query_scaffold_id_column_name,
                                                                   query_start_column_name,
                                                                   query_end_column_name,
                                                                   strand_column_name,
                                                                   args.inverted_scaffold_label)
    if genome_index > 0:  # apply listed inversion for all target genomes, i.e. for all except the first genome
        #print(genome_orderlist[genome_index - 1])
        #print(synteny_dict[genome_orderlist[genome_index - 1]])
        logging.info(TAB * 2 * TAB_NUM_PER_LVL + "Inverting target coordinates in synteny file...")
        synteny_dict[genome_orderlist[genome_index - 1]] = invert_coordinates_in_synteny_table(synteny_dict[genome_orderlist[genome_index - 1]],
                                                                                               config_dict[genome]["scaffold"]["invertlist"],
                                                                                               config_dict[genome]["scaffold"]["len_df"],
                                                                                               target_scaffold_id_column_name,
                                                                                               target_start_column_name,
                                                                                               target_end_column_name,
                                                                                               strand_column_name,
                                                                                               args.inverted_scaffold_label)

    for datatype in "scaffold", "superscaffold":
        if config_dict[genome][datatype]["centromere_df"] is not None:
            config_dict[genome][datatype]["centromere_df"] = invert_coordinates_in_region_table(config_dict[genome][datatype]["centromere_df"],
                                                                                                config_dict[genome][datatype]["invertlist"],
                                                                                                config_dict[genome][datatype]["len_df"],
                                                                                                "scaffold_id", "start", "end",
                                                                                                inverted_scaffolds_label=args.inverted_scaffold_label)
        if config_dict[genome][datatype]["invertlist"] is not None:
            invert_rename_dict = dict(zip(config_dict[genome][datatype]["invertlist"],
                                          [scaf + args.inverted_scaffold_label for scaf in config_dict[genome][datatype]["invertlist"]]))
            if config_dict[genome][datatype]["len_df"] is not None:
                config_dict[genome][datatype]["len_df"].rename(index=invert_rename_dict, inplace=True)
            if config_dict[genome][datatype]["orderlist"] is not None:
                config_dict[genome][datatype]["orderlist"].replace(invert_rename_dict, inplace=True)
            if config_dict[genome][datatype]["color_df"] is not None:
                config_dict[genome][datatype]["color_df"].rename(index=invert_rename_dict, inplace=True)
            config_dict[genome][datatype]["settings_df"].rename(index=invert_rename_dict, inplace=True)
#for genome in genome_orderlist:


#--------------------------------------------------------------------------------------

default_hex_color = mpl.colors.cnames["lightgrey"]
inversion_hex_color = mpl.colors.cnames[args.inversion_color] if args.inversion_color != "default" else default_hex_color
translocation_hex_color = mpl.colors.cnames[args.translocation_color] if args.translocation_color != "default" else default_hex_color

long_block_hex_color = "#00FF00"
short_block_hex_color = "#F4EA56"

for index, genome in zip(range(0, len(genome_orderlist) - 1), genome_orderlist[:-1]):
    target_genome = genome_orderlist[index + 1]
    datatype = "scaffold"
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
        logging.warn_scr(f"WARNING!!! Excel has a hardlimit of 31 char for sheetname. '{sheet_name}' is longer.Cutting to 31 chars...")
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
    for scaffold_id in config_dict[genome][datatype]["color_df"].index:
        query_scaffold_format_dict[scaffold_id] = workbook.add_format({'bg_color': config_dict[genome][datatype]["color_df"].loc[scaffold_id, "color"]})

    for scaffold_id in config_dict[target_genome][datatype]["color_df"].index:
        target_scaffold_format_dict[scaffold_id] = workbook.add_format({'bg_color': config_dict[target_genome][datatype]["color_df"].loc[scaffold_id, "color"]})

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
interchr_space_fraction = args.interscaf_space_fraction # 0.3

maximal_x = max_genome_length * (1 + interchr_space_fraction)

height = args.chromosome_height
length = 1000000
distance = args.genome_distance

maximal_y = height * len(genome_orderlist) + distance * (len(genome_orderlist) - 1)

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


def add_chromosome_patch_function(row, genome, datatype): #centromere_df):
    #print(row)
    #print("AAAAA")
    #print(pd.isna(row.iloc[5]))
    #print("BBBBBBBBbb")
    #print((row.iloc[5] + row.iloc[1]) if not pd.isna(row.iloc[5]) else None)
    #print((row.iloc[6] + row.iloc[1]) if not pd.isna(row.iloc[5]) else None)
    #print(row.name)
    return LinearChromosome(row.iloc[1], row.iloc[2], row.iloc[0], height,
                            rounded=config_dict[genome][datatype]["settings_df"].loc[row.name, "rounded"],
                            stranded=config_dict[genome][datatype]["settings_df"].loc[row.name, "stranded"],
                            stranded_end=config_dict[genome][datatype]["settings_df"].loc[row.name, "stranded_end"],
                            x_scale_factor=args.smooth_multiplicator * maximal_x/10 / maximal_y, #maximal_x/10 / maximal_y,
                            zorder=zorder_dict["chromosomes"],
                            edgecolor=config_dict[genome][datatype]["settings_df"].loc[row.name, "edge_color"],
                            facecolor=config_dict[genome][datatype]["settings_df"].loc[row.name, "fill_color"],
                            alpha=0.9,
                            linewidth=0.3,
                            centromere_start=(row.iloc[5]) if not pd.isna(row.iloc[5]) else None, #+ row.iloc[1]
                            centromere_end=(row.iloc[6]) if not pd.isna(row.iloc[5]) else None, #+ row.iloc[1]
                            show_centromere=True)


def add_chromosome_line_function(row, height):
    return Line2D(xdata=(row.iloc[1], row.iloc[1] + row.iloc[0]),
                  ydata=(row.iloc[2] + height/2, row.iloc[2] + height/2,),
                  color="black",
                  zorder=zorder_dict["chromosome_lines"],
                  alpha=0.4,
                  linewidth=0.3,)




for genome, index, genome_color, genome_label in zip(genome_orderlist, range(0, len(genome_orderlist)), color_list, genome_labellist): #genome_orderlist):
    #print(centromere_df_dict[genome])
    datatype = "scaffold"
    interchr_space = ((maximal_x - total_len_dict[genome]) / (chr_number_dict[genome] - 1)) if chr_number_dict[genome] > 1 else 0
    config_dict[genome][datatype]["len_df"]["x_offset"] = (config_dict[genome][datatype]["len_df"]["length"].cumsum().shift(periods=1, fill_value=0) + np.array(range(0, chr_number_dict[genome])) * interchr_space)
    config_dict[genome][datatype]["len_df"]["y_offset"] = (height + distance) * index
    config_dict[genome][datatype]["len_df"]["color"] = default_genome_color_dict[genome]

    config_dict[genome][datatype]["len_df"]["label"] = pd.Series(list(config_dict[genome][datatype]["len_df"].index),
                                                                 index=config_dict[genome][datatype]["len_df"].index).apply(lambda s: s[args.scaffold_prefix_cut:])
    config_dict[genome][datatype]["len_df"]["centromere_start"] = config_dict[genome][datatype]["centromere_df"]["start"]
    config_dict[genome][datatype]["len_df"]["centromere_end"] = config_dict[genome][datatype]["centromere_df"]["end"]
    #print(config_dict[genome]["scaffold"]["len_df"])
    patch_collection = PatchCollection(config_dict[genome]["scaffold"]["len_df"].apply(partial(add_chromosome_patch_function,
                                                                                               genome=genome,
                                                                                               datatype="scaffold"),
                                                                                       axis=1),
                                       match_original=True,
                                       antialiased=False,
                                       zorder=zorder_dict["chromosomes"])
    ax.add_collection(patch_collection)

    for line in config_dict[genome]["scaffold"]["len_df"].apply(partial(add_chromosome_line_function, height=height), axis=1):
        ax.add_line(line)

    if not args.hide_chromosome_labels:
        for scaffold in config_dict[genome]["scaffold"]["len_df"].index:
            #print(scaffold)
            #print(config_dict[genome]["scaffold"]["settings_df"])
            if config_dict[genome]["scaffold"]["settings_df"].loc[scaffold, "show_label"]:
                ax.annotate(config_dict[genome]["scaffold"]["len_df"].loc[scaffold, "label"],
                            xy=(config_dict[genome]["scaffold"]["len_df"].loc[scaffold, "x_offset"] + config_dict[genome]["scaffold"]["len_df"].loc[scaffold, "length"]/2,
                                config_dict[genome]["scaffold"]["len_df"].loc[scaffold, "y_offset"] + height), xycoords='data',
                            fontsize=config_dict[genome]["scaffold"]["settings_df"].loc[scaffold, "label_fontsize"],
                            xytext=(0, 0), textcoords='offset points',
                            rotation=config_dict[genome]["scaffold"]["settings_df"].loc[scaffold, "label_angle"],
                            ha="center", va="bottom",
                            color="black",
                            zorder=zorder_dict["label"])

    if config_dict[genome]["genome"]["show_label"]:
        ax.annotate(genome_label,
                    xy=(- border_offset_fraction / 2 * maximal_x, (height + distance) * index + height/2),
                    xycoords='data',
                    fontsize=config_dict[genome]["genome"]["label_fontsize"],  #args.genome_label_fontsize,
                    fontstyle="italic",
                    xytext=(0, 0), textcoords='offset points', rotation=config_dict[genome]["genome"]["label_angle"],
                    ha="right", va="center",
                    color="black",
                    zorder=zorder_dict["label"])


def connector_function(row, length_df_dict, top_genome, bottom_genome, default_color="lightgrey", connector_color_idx=8,
                      top_scaffold_idx=3, top_start_idx=4, top_end_idx=5, bottom_scaffold_idx=0, bottom_start_idx=1, bottom_end_idx=2, strand_idx=6):
    con_len = (connector_color_idx + 1) if connector_color_idx is not None else None
    y_chr_shift = height / 2
    #print("AAAAAAAAA")
    #print(top_genome)
    #print(length_df_dict[top_genome])
    #print("BBBBBBBBB")
    #print(bottom_genome)
    #print(length_df_dict[bottom_genome])
    return CubicBezierConnector(
                                 (row.iloc[top_start_idx] + length_df_dict[top_genome].loc[row.iloc[top_scaffold_idx], "x_offset"], length_df_dict[top_genome].loc[row.iloc[top_scaffold_idx], "y_offset"] + y_chr_shift),
                                 (row.iloc[top_end_idx] + length_df_dict[top_genome].loc[row.iloc[top_scaffold_idx], "x_offset"], length_df_dict[top_genome].loc[row.iloc[top_scaffold_idx], "y_offset"] + y_chr_shift),

                                ((row.iloc[bottom_start_idx] if row.iloc[strand_idx] == "+" else row.iloc[bottom_end_idx]) + length_df_dict[bottom_genome].loc[row.iloc[bottom_scaffold_idx], "x_offset"],
                                 length_df_dict[bottom_genome].loc[row.iloc[bottom_scaffold_idx], "y_offset"] + y_chr_shift),

                                ((row.iloc[bottom_end_idx] if row.iloc[strand_idx] == "+" else row.iloc[bottom_start_idx]) + length_df_dict[bottom_genome].loc[row.iloc[bottom_scaffold_idx], "x_offset"],
                                 length_df_dict[bottom_genome].loc[row.iloc[bottom_scaffold_idx], "y_offset"] + y_chr_shift),

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

for genome_index in range(0, len(genome_orderlist) - 1):
    query_genome = genome_orderlist[genome_index]
    target_genome = genome_orderlist[genome_index + 1]
    len_df_dict = {
                   query_genome: config_dict[query_genome]["scaffold"]["len_df"],
                   target_genome: config_dict[target_genome]["scaffold"]["len_df"]
                   }
    if "connector_zorder" in synteny_dict[query_genome]:
        synteny_dict[query_genome]["connector_zorder"] += zorder_dict["connector"]
        connector_collection_dict[query_genome] = {}
        for zorder in sorted(synteny_dict[query_genome]["connector_zorder"].unique()):
            #print(config_dict[genome]["scaffold"]["len_df"])
            connector_collection_dict[query_genome][zorder] = PatchCollection(synteny_dict[query_genome][synteny_dict[query_genome]["connector_zorder"] == zorder].apply(partial(connector_function,
                                                                                                                                                                       length_df_dict=len_df_dict,
                                                                                                                                                                       top_genome=target_genome,
                                                                                                                                                                       connector_color_idx=connector_color_idx,
                                                                                                                                                                       bottom_genome=query_genome,
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
            ax.add_collection(connector_collection_dict[query_genome][zorder])
    else:
        connector_collection_dict[query_genome] = PatchCollection(synteny_dict[query_genome].apply(partial(connector_function,
                                                                                                           length_df_dict=len_df_dict,
                                                                                                           top_genome=target_genome,
                                                                                                           connector_color_idx=connector_color_idx,
                                                                                                           bottom_genome=genome,
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
        ax.add_collection(connector_collection_dict[query_genome])

plt.subplots_adjust(left=args.subplots_adjust_left, right=args.subplots_adjust_right, bottom=args.subplots_adjust_bottom,
                    top=args.subplots_adjust_top)
plt.title(args.title, fontsize=args.title_fontsize)
for ext in args.output_formats:
    if args.manual_figure_adjustment:
        plt.savefig(f"{args.output_prefix}.{ext}")
    else:
        plt.savefig(f"{args.output_prefix}.{ext}", bbox_inches="tight",)
