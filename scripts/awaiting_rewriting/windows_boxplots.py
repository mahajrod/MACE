#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from _collections import OrderedDict


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
                         "Could be not set. "
                         "Each X-chromosome could contain coordinates of "
                         "pseudoautosomal region (PAR) following this scheme: <scaffold>:<PARstart>-<PARend>."
                         "PAR coordinates must be 1-based")

parser.add_argument("-w", "--window_size", action="store", dest="window_size", required=True, type=int,
                    help="Size of the windows. Required.")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", required=True, type=int,
                    help="Step of the windows. Required.")
parser.add_argument("-m", "--multiplier", action="store", dest="multiplier", default=1000, type=int,
                    help="Multiplier for counts. Default:1000, ie counts will scaled per 1000 bp")
parser.add_argument("-y", "--ylabel", action="store", dest="ylabel", default="variants/kbp",
                    help="Label for Y axis. Default: 'variants/kbp'")
parser.add_argument("--ymax", action="store", dest="ymax", type=float,
                    help="Max value to show for Y axis. Not set")
"""

parser.add_argument("--figure_height_per_scaffold", action="store", dest="figure_height_per_scaffold",
                    type=float, default=0.5,
                    help="Height of figure per chromosome track. Default: 0.5")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Print additional info to stdout")
"""

parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float, default=0.15,
                    help="Adjust left border of subplots on the figure. Default: 0.17")
parser.add_argument("--subplots_adjust_top", action="store", dest="subplots_adjust_top", type=float, default=0.98,
                    help="Adjust top border of subplots on the figure. Default: 0.98")
parser.add_argument("--subplots_adjust_right", action="store", dest="subplots_adjust_right", type=float, default=0.98,
                    help="Adjust right border of subplots on the figure. Default: 0.98")
parser.add_argument("--subplots_adjust_bottom", action="store", dest="subplots_adjust_bottom", type=float, default=0.25,
                    help="Adjust bottom border of subplots on the figure. Default: 0.17")
args = parser.parse_args()


def parse_x_string(x_list):
    tmp = list(map(lambda s: s.split(":"), x_list))
    parsed_dict_list = []
    for element in tmp:
        element_dict = OrderedDict({"id": element[0]})
        if len(element) > 1:
            element_dict["PARstart"], element_dict["PARend"] = list(map(int, element[1].split("-")))
            element_dict["PARstart"] -= 1
        parsed_dict_list.append(element_dict)

    return parsed_dict_list

number_of_files = len(args.input_list)

if number_of_files != len(args.label_list):
    raise ValueError("ERROR!!! Number of input files({0}) is not equal to number of labels({1})".format(number_of_files, len(args.label_list)))

if args.x_list is None:
    x_chr_list = []
elif len(args.x_list) == 1:
    x_chr_list = parse_x_string(args.x_list) * number_of_files
else:
    if number_of_files != len(args.x_list):
        raise ValueError("ERROR!!! Number of X chromosome ids is required to be one or equal to number of input files({0})".format(number_of_files))
    else:
        x_chr_list = parse_x_string(args.x_list)

file_dict = OrderedDict(zip(args.label_list, args.input_list))
x_chr_dict = OrderedDict(zip(args.label_list, x_chr_list)) if x_chr_list else OrderedDict()

df_dict = OrderedDict({"all": OrderedDict(),
                       "noX": OrderedDict(),
                       "onlyX": OrderedDict(),
                       "noPAR": OrderedDict(),
                       "PAR": OrderedDict()})

df_number_dict = OrderedDict()

for entry in file_dict:
    df_dict["all"][entry] = pd.read_csv(file_dict[entry], sep="\t", index_col=["CHROM",])
    df_dict["all"][entry]["density"] = df_dict["all"][entry]["All"] * args.multiplier / args.window_size
    df_number_dict[entry] = 1

    if x_chr_dict:
        df_dict["noX"][entry] = df_dict["all"][entry].loc[df_dict["all"][entry].index != x_chr_dict[entry]["id"], :]
        df_dict["onlyX"][entry] = df_dict["all"][entry].loc[df_dict["all"][entry].index == x_chr_dict[entry]["id"], :]
        df_number_dict[entry] += 2
        if "PARstart" in x_chr_dict[entry]:
            PAR_boolean_windows = (df_dict["onlyX"][entry]["WINDOW"] * args.window_step >= x_chr_dict[entry]["PARstart"]) & (
                        df_dict["onlyX"][entry]["WINDOW"] * args.window_step + args.window_size <= x_chr_dict[entry]["PARend"])
            df_dict["noPAR"][entry] = df_dict["onlyX"][entry][~PAR_boolean_windows]
            df_dict["PAR"][entry] = df_dict["onlyX"][entry][PAR_boolean_windows]
            df_number_dict[entry] += 2
    #if x_chr_dict:
    #    df_dict["noX"][entry]["density"] = df_dict["noX"][entry]["All"] * args.multiplier / args.window_size
    #    df_dict["onlyX"][entry]["density"] = df_dict["onlyX"][entry]["All"] * args.multiplier / args.window_size

df_number_list = list(df_number_dict.values())

figure_height = 4
figure_width = max(1, int(sum(df_number_list) / 2))
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
        for label in df_dict:
            if not df_dict[label]:
                continue
            data_list.append(df_dict[label][entry]["density"])
            label_list.append(label)
            if not x_pos_list:
                x_pos_list.append(0)
            else:
                x_pos_list.append(x_pos_list[-1] + distance)
                distance = 1
        distance = otter_distance
    plt.xticks(rotation=45)
else:
    distance = 1
    for entry in file_dict:
        data_list.append(df_dict["all"][entry]["density"])
        label_list.append(entry)
        if not x_pos_list:
            x_pos_list.append(0)
        else:
            x_pos_list.append(x_pos_list[-1] + distance)
    plt.xticks(rotation=45, fontstyle='italic')

plt.boxplot(data_list, labels=label_list, positions=x_pos_list)
plt.ylabel(args.ylabel)
plt.ylim(ymin=-0.1)
plt.grid(linestyle='dotted')

max_data_y = max(list(map(max, data_list)))
min_data_y = min(list(map(min, data_list)))
plt.ylim(ymax=args.ymax)
"""
ticks = [tick for tick in plt.gca().get_xticklabels()]
for i, t in enumerate(ticks):
    print(t.get_horizontalalignment())
    print(t.get_va())
    print(t.get_position())
    print(t.get_bbox_patch())
    print("Label ", i, ", data: ", t.get_text(), " ; ", t.get_window_extent())
"""
if x_chr_dict:
    prev = 0
    for entry in df_number_dict:
        if x_chr_dict:
            label_x_coord = (x_pos_list[prev] + x_pos_list[prev + df_number_dict[entry] - 1])/2
        else:
            label_x_coord = x_pos_list[prev]
        plt.text(label_x_coord, min_data_y - (args.ymax if args.ymax else max_data_y)/3, entry,
                 fontstyle="italic", horizontalalignment='center')
        prev += df_number_dict[entry]

plt.subplots_adjust(top=args.subplots_adjust_top, bottom=args.subplots_adjust_bottom,
                    left=args.subplots_adjust_left, right=args.subplots_adjust_right)
for ext in "png", "svg":
    plt.savefig("{0}.{1}".format(args.output_prefix, ext), transparent=False)
