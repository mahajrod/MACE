#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from collections import OrderedDict
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input tab file with at least two columns containing x and y values."
                         "Default: stdin")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Column separator in the input file"
                         "Default: '\t'")
parser.add_argument("-r", "--header_row", action="store", dest="header_row", default=None,
                    help="Index(zero-based) of header row in the input file."
                         "Default: header is absent")
parser.add_argument("-c", "--comments_prefix", action="store", dest="comments_prefix", default=None,
                    help="Prefix of the comment rows in the input file."
                         "Default: not set")
parser.add_argument("-x", "--x_column_index", action="store", dest="x_column_index", default=0, type=int,
                    help="Index (zero-based) of column in the input file containing x values."
                         "Default: 0")
parser.add_argument("-y", "--y_column_index", action="store", dest="y_column_index", default=1, type=int,
                    help="Index (zero-based) of column in the input file containing y values."
                         "Default: 1")
parser.add_argument("-w", "--weights_column_index", action="store", dest="weights_column_index", default=None, type=int,
                    help="Index (zero-based) of column in the input file containing weights."
                         "Default: absent")
parser.add_argument("--x_label", action="store", dest="x_label",
                    help="Label for x axis. Default: not set")
parser.add_argument("--y_label", action="store", dest="y_label",
                    help="Label for y axis. Default: not set")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=300,
                    help="Dpi of figure")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=("png", ),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg,png")

parser.add_argument("-l", "--title", action="store", dest="title", default="Variant density",
                    help="Suptitle of figure. Default: 'Variant density'")

parser.add_argument("--ymin", action="store", dest="ymin", type=float, default=None,
                    help="Minimum limit for Y axis . Default: -0.1")
parser.add_argument("--ymax",  action="store", dest="ymax", type=float, default=None,
                    help="Maximum limit for Y axis. Default: not set")

parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float, default=0.05,
                    help="Adjust left border of subplots on the figure. Default: 0.05")
parser.add_argument("--subplots_adjust_top", action="store", dest="subplots_adjust_top", type=float, default=0.95,
                    help="Adjust top border of subplots on the figure. Default: 0.05")
parser.add_argument("--subplots_adjust_right", action="store", dest="subplots_adjust_right", type=float, default=0.95,
                    help="Adjust right border of subplots on the figure. Default: 0.95")
parser.add_argument("--subplots_adjust_bottom", action="store", dest="subplots_adjust_bottom", type=float, default=0.05,
                    help="Adjust bottom border of subplots on the figure. Default: 0.95")
parser.add_argument("--subplots_adjust_hspace", action="store", dest="subplots_adjust_hspace", type=float, default=0.1,
                    help="Adjust height space between subplots on the figure. Default: 0.1")
parser.add_argument("--subplots_adjust_wspace", action="store", dest="subplots_adjust_wspace", type=float, default=0.1,
                    help="Adjust width space between  subplots on the figure. Default: 0.1")
args = parser.parse_args()

data_df = pd.read_csv(args.input, sep=args.separator,
                      usecols=[args.x_column_index, args.y_column_index, args.weights_column_index] if args.weights_column_index else [args.x_column_index, args.y_column_index],
                      names=["x", "y", "w"] if args.weights_column_index else ["x", "y"],
                      comment=args.comments_prefix, header=args.header_row)

fig, ax_array = plt.subplots(2, 2, dpi=args.dpi, figsize=(12, 12), width_ratios=(2, 1), height_ratios=(1, 2))
ax_array[1][0].sharex(ax_array[0][0])
ax_array[1][0].sharey(ax_array[1][1])

x_bins = []
y_bins = []

ax_array[0][0].hist(data_df["x"], weights=data_df["w"] if args.weights_column_index else None,
                    bins=x_bins, )
ax_array[0][0].set_ylabel(args.x_label if args.y_label else "Y")

ax_array[1][0].hist2d(data_df["x"], data_df["y"], weights=data_df["w"] if args.weights_column_index else None,
                      bins=[x_bins, y_bins],
                      norm=mpl.colors.LogNorm() if color_scale == "log" else None)

ax_array[1][0].set_xlabel(args.x_label if args.x_label else "X")
ax_array[1][0].set_ylabel(args.y_label if args.y_label else "Y")

ax_array[1][1].hist(data_df["y"], weights=data_df["w"] if args.weights_column_index else None,
                    bins=y_bins, orientation="horizontal")
ax_array[1][0].set_xlabel(args.x_label if args.x_label else "X")

plt.subplots_adjust(left=args.subplots_adjust_left, right=args.subplots_adjust_right,
                    top=args.subplots_adjust_top, bottom=args.subplots_adjust_bottom,
                    hspace=args.subplots_adjust_hspace, wspace=args.subplots_adjust_wspace)





