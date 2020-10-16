#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with two columns containing label in the first one and filename in the second."
                         "Boxplots will be drawn in the same order as labels")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=300,
                    help="Dpi of figure")

parser.add_argument("-f", "--size_of_figure", action="store", dest="size_of_figure", type=lambda s: s.split(","),
                    default=(6, 12),
                    help="Size of figure in inches. X and Y values should be separated by comma. Default: 6,12")

parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=("png", ),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg,png")

parser.add_argument("-l", "--title", action="store", dest="title", default="Variant density",
                    help="Suptitle of figure. Default: 'Variant density'")

parser.add_argument("--denominator", action="store", dest="denominator", default=1000, type=float,
                    help="Denominator for variant counts. "
                         "Default: 1000, i.e variant counts will be scaled to per 1 kbp ")

parser.add_argument("--ymin", action="store", dest="ymin", type=float, default=-0.1,
                    help="Minimum limit for Y axis . Default: -0.1")
parser.add_argument("--ymax",  action="store", dest="ymax", type=float, default=None,
                    help="Maximum limit for Y axis. Default: not set")
parser.add_argument("--yticklist",  action="store", dest="yticklist", type=lambda s: list(map(float, s.split())),
                    default=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2, 3, 4, 5],
                    help="Comma-separated tick list for Y axis. "
                         "Default: 0.05,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2,3,4,5")
parser.add_argument("--rotation",  action="store", dest="rotation", type=float, default=90,
                    help="Rotation angle for X labels. Default: 90")

parser.add_argument("--horizontal_lines",  action="store", dest="horizontal_lines",
                    type=lambda s: list(map(float, s.split())),
                    help="Comma-separated list of y-coordinates to draw horizontal lines. "
                         "Default: not set")
"""
parser.add_argument("-q", "--figure_width", action="store", dest="figure_width", default=12, type=int,
                    help="Width of figure in inches. Default: 12")
parser.add_argument("-u", "--figure_height_scale_factor", action="store", dest="figure_height_scale_factor",
                    default=0.5, type=float,
                    help="Figure height scale factor. Figure height is calculated in inches as "
                         "int(figure_scale_factor * scaffold_number * sample_number). Default: 0.5")
"""


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

args = parser.parse_args()

with open(args.input, "r") as in_fd:
    file_dict = OrderedDict([line.strip().split("\t") for line in in_fd ])

df_dict = OrderedDict({})
for entry in file_dict:
    df_dict[entry] = pd.read_csv(file_dict[entry], sep="\t", index_col=["CHROM",])
    df_dict[entry]["density"] = df_dict[entry]["All"] / args.denominator

fig, ax = plt.subplots(figsize=args.size_of_figure, dpi=args.dpi)
plt.xticks(rotation=args.rotation)
plt.boxplot([df_dict[entry]["density"] for entry in df_dict], labels=list(df_dict.keys()))


plt.yticks(args.yticklist)
if args.horizontal_lines:
    for ycoord in args.horizontal_lines:
        plt.axhline(y=ycoord, color="red", linestyle="--", linewidth=0.5)
plt.ylim(ymax=args.ymax, ymin=args.ymin)
plt.subplots_adjust(left=args.subplots_adjust_left, right=args.subplots_adjust_right,
                    top=args.subplots_adjust_top, bottom=args.subplots_adjust_bottom)

plt.title(args.title)
for ext in args.output_formats:
    plt.savefig("{0}.{1}".format(args.output_prefix, ext))



