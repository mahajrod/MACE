#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from _collections import OrderedDict

from RouToolPa.Collections.General import SynDict, IdList
from MACE.Routines import Visualization, StatsVCF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_list", action="store", dest="input_list", required=True, type=lambda s: s.split(","),
                    help="Comma-separated list of files with variant counts in windows.")
parser.add_argument("-l", "--label_list", action="store", dest="label_list", required=True,type=lambda s: s.split(","),
                    help="Comma-separated list of labels for files. Should has the same length as input_list")

parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=("png", "svg"),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg, png")

parser.add_argument("-x", "--x_list", action="store", dest="x_list", type=lambda s: s.split(","),
                    default=None,
                    help="Comma-separated list of ids of X chromosomes. "
                         "Should include either one element if X chromosome id is thesame for all input files or has the same length as input list."
                         "Could be not set ")

parser.add_argument("-w", "--window_size", action="store", dest="window_size", required=True, type=int,
                    help="Size of the windows. Required.")
parser.add_argument("-m", "--multiplier", action="store", dest="multiplier", default=1000, type=int,
                    help="Multiplier for counts. Default:1000, ie counts will scaled per 1000 bp")
parser.add_argument("-y", "--ylabel", action="store", dest="ylabel", default="variants/kbp",
                    help="Label for Y axis. Default: 'variants/kbp'")
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
                    type=lambda s: IdList(filename=s) if os.path.exists(s) else s.split(","),
                    help="Comma-separated list of the only scaffolds to draw. Default: all")

parser.add_argument("-b", "--scaffold_black_list", action="store", dest="scaffold_black_list", default=[],
                    type=lambda s: IdList(filename=s) if os.path.exists(s) else s.split(","),
                    help="Comma-separated list of scaffolds to skip at drawing. Default: not set")

parser.add_argument("-y", "--sort_scaffolds", action="store_true", dest="sort_scaffolds", default=False,
                    help="Order  scaffolds according to their names. Default: False")

parser.add_argument("-z", "--scaffold_ordered_list", action="store", dest="scaffold_ordered_list", default=[],
                    type=lambda s: IdList(filename=s) if os.path.exists(s) else s.split(","),
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
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Print additional info to stdout")
"""
args = parser.parse_args()

number_of_files = len(args.input_list)

if number_of_files != len(args.label_list):
    raise ValueError("ERROR!!! Number of input files({0}) is not equal to number of labels({1})".format(number_of_files, len(args.label_list)))

if args.x_list is None:
    x_chr_list = []
elif len(args.x_list) == 1:
    x_chr_list = args.x_list * number_of_files
else:
    if number_of_files != len(args.x_list):
        raise ValueError("ERROR!!! Number of X chromosome ids is required to be one or equal to number of input files({0})".format(number_of_files))
    else:
        x_chr_list = args.x_list

file_dict = OrderedDict(zip(args.label_list, args.input_list))
x_chr_dict = OrderedDict(zip(args.label_list, x_chr_list)) if x_chr_list else OrderedDict()

df_dict = OrderedDict({})
noX_df_dict = OrderedDict({})
onlyX_df_dict = OrderedDict({})

for entry in file_dict:
    df_dict[entry] = pd.read_csv(file_dict[entry], sep="\t", index_col=["CHROM",])
    if x_chr_dict:
        noX_df_dict[entry] = df_dict[entry].loc[df_dict[entry].index != x_chr_dict[entry], :]
        onlyX_df_dict[entry] = df_dict[entry].loc[df_dict[entry].index == x_chr_dict[entry], :]

    df_dict[entry]["density"] = df_dict[entry]["All"] * args.multiplier / args.window_size
    if x_chr_dict:
        noX_df_dict[entry]["density"] = noX_df_dict[entry]["All"] * args.multiplier / args.window_size
        onlyX_df_dict[entry]["density"] = onlyX_df_dict[entry]["All"] * args.multiplier / args.window_size

figure_height = 6
figure_width = (number_of_files * 3) if x_chr_dict else number_of_files
dpi = 300
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figure_width, figure_height), dpi=dpi)
fig.patch.set_color('white')

data_list = []
label_list = []
x_pos_list = []

inner_distance = 1
otter_distance = 2


if x_chr_dict:
    distance = inner_distance
    for entry in file_dict:
        for data, label in zip((df_dict, noX_df_dict, onlyX_df_dict), ("all", "noX", "onlyX")):
            data_list.append(data[entry]["density"])
            label_list.append(label)
            if not x_pos_list:
                x_pos_list.append(0)
            else:
                x_pos_list.append(x_pos_list[-1] + distance)
                distance = 1
        distance = otter_distance
    plt.xticks(rotation=0)
else:
    distance = 1
    for entry in file_dict:
        data_list.append(df_dict[entry]["density"])
        label_list.append(entry)
        if not x_pos_list:
            x_pos_list.append(0)
        else:
            x_pos_list.append(x_pos_list[-1] + distance)
    plt.xticks(rotation=45)


plt.boxplot(data_list, labels=label_list, positions=x_pos_list)
plt.ylabel(args.ylabel)
plt.ylim(ymin=-0.1)
plt.grid(linestyle='dotted')


max_data_y = max(list(map(max, data_list)))
min_data_y = min(list(map(min, data_list)))
if x_chr_dict:
    for index in range(0, number_of_files):
        if x_chr_dict:
            label_x_coord = x_pos_list[index * number_of_files] + (x_pos_list[(index + 1) * number_of_files - 1] - x_pos_list[index * number_of_files])/2
        else:
            label_x_coord = x_pos_list[index]
        plt.text(label_x_coord, min_data_y - max_data_y/6 , args.label_list[index],
                 fontstyle="italic", horizontalalignment='center')

plt.subplots_adjust(top=0.98, bottom=0.17)
for ext in "png", "svg":
    plt.savefig("{0}.{1}".format(args.output_prefix, ext), transparent=False)
