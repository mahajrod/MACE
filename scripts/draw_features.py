#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

import pandas as pd

from RouToolPa.Parsers.STR import CollectionSTR
from RouToolPa.Parsers.GFF import CollectionGFF
from RouToolPa.Parsers.BLAST import CollectionBLAST
from RouToolPa.Parsers.BED import CollectionBED

from MACE.Routines import Visualization
from MACE.Routines import Parsing

from MACE.Routines.Parsing import NewlinePreservingArgParserHelpFormatter


parser = argparse.ArgumentParser(formatter_class=NewlinePreservingArgParserHelpFormatter)
# ---- Input/Output options ----

# -------- Main Input options --------
parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with selected features")
parser.add_argument("-t", "--input_type", action="store", dest="input_type", default="bed",
                    help="Type of input file. Allowed: 'bed' (default), 'bed_colored', 'bed_track', 'bedgraph', 'bedgraph_colored', "
                         " 'bed_with_header', 'table', 'table_colored', 'gff' 'tab6', 'tab6_colored', 'STR'. All bed-like formats must have no track lines.\n"
                         "\t'bed': BED format. All columns except first three are ignored. Comment lines are allowed and must start from '#'. NO HEADER!\n"
                         "\t'bed_colored': BED format with four columns: 'scaffold', 'start', 'end', 'color'. Comment lines are not allowed. NO HEADER!\n"
                         "\t'bed_track': BED format with five columns: 'scaffold', 'start', 'end', 'value', 'color'. Comment lines are not allowed. NO HEADER!\n"
                         "\t'bedgraph': a synonym to 'bed'. Fourth column of bedgraph is ignored.\n"
                         "\t'bedgraph_colored': a synonym to 'bed_track'.\n"
                         "\t'bed_with_header': 'BED' format without track lines, but WITH a header line. Not a true BED!\n"
                         "\t'table': a tab separated table WITH a header line. Scaffold, start and end columns must be set via script options. "
                            "Other columns are ignored."
                            "Comment lines are allowed and must start from '#'.\n"
                            "If it is one-based like GFF, set '--one_based_coordinates' script option. Otherwise it will be treated as zero-based like BED."
                         "\t'table_colored': A tab separated table with a header line. Scaffold, start, end and color columns must be set via script options. "
                            "Other columns are ignored."
                            "Comment lines are  NOT ALLOWED.\n"
                            "If it is one-based like GFF, set '--one_based_coordinates' script option. Otherwise it will be treated as zero-based like BED."
                         "\t'gff': GFF format. Only coordinates are used. Comment lines are allowed and must start from '#'.\n"
                         "\t'gtf': GTF format. Only coordinates are used. Comment lines are allowed and must start from '#'.\n"
                         "\t'blast6': BLAST output 6 format.\n"
                         "\t'blast6_colored': BLAST output 6 format with header and additional color column. Set color column via script options. Names of other columns are ignored\n"
                         "\t'STR': a very specific format related to STRs and in silico PCR. Highly likely you don't need it")

parser.add_argument("--scaffold_column_name", action="store", dest="scaffold_column_name", default="scaffold",
                    help="Name of column in feature file with scaffold ids . Default: dependent on format of the file")
parser.add_argument("--start_column_name", action="store", dest="start_column_name", default="start",
                    help="Name of column in feature file with starts. Default: dependent on format of the file")
parser.add_argument("--end_column_name", action="store", dest="end_column_name", default="end",
                    help="Name of column in feature file with ends. Default: dependent on format of the file")
parser.add_argument("--color_column_name", action="store", dest="color_column_name",  default="color",
                    help="Name of column in feature file with color. Default: not set")
parser.add_argument("--one_based_coordinates", action="store_true", dest="one_based_coordinates", default=False,
                    help="Input uses one-based coordinate system (like GFF format). "
                         "This option affects only 'table' input format and its derivatives like 'table_colored'. "
                         "For other formats input coordinate system is set according to the specification of corresponding format. "
                         "Internally the script and related libraries use zero-based coordinate system (like in BED format), which is also set as a default. "
                         "Default: False")

parser.add_argument("-n", "--scaffold_length_file", action="store", dest="scaffold_length_file", required=True,
                    help="File with lengths of scaffolds")

parser.add_argument("--centromere_bed", action="store", dest="centromere_bed", required=False,
                    type=str, help="BED file with coordinates of centromeres. Optional.")
# -------- End of Main Input options --------

# -------- Input Filtering, Renaming, Sorting and Highlighting options --------
parser.add_argument("-a", "--scaffold_whitelist", action="store", dest="scaffold_whitelist",
                    help="Comma-separated list of the only scaffolds to draw. Default: not set, number of longest scaffolds to show is controlled by '--max_scaffolds'")
parser.add_argument("--max_scaffolds", action="store", dest="max_scaffolds",
                    default=50,
                    type=int,
                    help="Maximal number of longest scaffolds from input file to show. This option works only if --scaffold_whitelist is not set. Default: 50")

parser.add_argument("-z", "--scaffold_orderlist", action="store", dest="scaffold_orderlist",
                    help="Comma-separated list of scaffolds to draw first and exactly in same order. "
                         "Scaffolds absent in this list are drawn last and in order according to vcf file . "
                         "Default: not set")

parser.add_argument("--scaffold_syn_file", action="store", dest="scaffold_syn_file",
                    help="File with scaffold id synonyms")
parser.add_argument("--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in synonym file. Default: 0")
parser.add_argument("--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in synonym file synonym. Default: 1")
parser.add_argument("--highlight_file", action="store", dest="highlight_file",
                    type=lambda s: pd.read_csv(s, header=0, index_col=0, sep="\t"),
                    help="Tab-separated file with two columns ('scaffold' and 'color'). "
                         "Scaffold ids are ids after renaming. Must contain header.")
# -------- End of Input Filtering, Renaming and Sorting options --------

# -------- Legend options --------
parser.add_argument("-g", "--legend", action="store", dest="legend",
                    help="File with legend for feature colors containing two columns with color and legend text")
parser.add_argument("--legend_colormap", action="store", dest="legend_colormap", default="jet",
                    help="Matplotlib colormap to use for legend. Default: jet")
# -------- End of Legend options --------

# -------- Output options -------
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=("png", "svg"),
                    help="Comma-separated list of formats (supported by matplotlib) for output figure. Default: svg,png")
# -------- End of Output options -------

# ---- End of Input/Output options ----

# ---- Drawing options ----

# -------- Feature options --------
parser.add_argument("--feature_shape", action="store", dest="feature_shape", default="rectangle",
                    help="Shape of features. Allowed: rectangle(default), circle, ellipse")
parser.add_argument("--default_color", action="store", dest="default_color", default="red",
                    help="Default color used for all features if color column is not set. Default: red")
parser.add_argument("--enforce_default_color", action="store_true", dest="enforce_default_color", default=False,
                    help="Enforce default color for all features even if input format natively has a color column. Default: False")
# -------- End of Feature options --------

# -------- Chromosome track options --------
parser.add_argument("--stranded", action="store_true", dest="stranded", default=False,
                    help="Stranded features and tracks. Default: False")
parser.add_argument("--rounded", action="store_true", dest="rounded", default=False,
                    help="Rounded tracks. Default: False")
parser.add_argument("--stranded_end", action="store_true", dest="stranded_end", default=False,
                    help="Stranded ends for tracks. Works only if --stranded is set. Default: False")
parser.add_argument("--fill_empty_tracks", action="store_true", dest="fill_empty_tracks", default=False,
                    help="Fill empty tracks (without features) with color set by --empty_color. Default: False")
parser.add_argument("--empty_color", action="store", dest="empty_color", default="lightgrey",
                    help="Color used to fill empty tracks. Ignored if --fill_empty_tracks is not set. "
                         "Default: 'lightgrey'")
# -------- End of Chromosome track options --------

# -------- Title options --------
parser.add_argument("-l", "--title", action="store", dest="title", default="Coverage",
                    help="Suptitle of figure. Default: 'Coverage'")
parser.add_argument("--title_fontsize", action="store", dest="title_fontsize", default=20, type=int,
                    help="Fontsize of the figure. Default: 20")
# -------- End of Title options --------

# -------- Label and Tick options --------
parser.add_argument("--hide_track_label", action="store_true", dest="hide_track_label", default=False,
                    help="Hide track label. Default: False")
parser.add_argument("--x_tick_type", action="store", dest="x_tick_type", default="nucleotide",
                    help="Type of xticks. Allowed: 'nucleotide' (default), 'int_number', 'float_number'")
parser.add_argument("--x_tick_fontsize", action="store", dest="x_tick_fontsize", type=int, default=None,
                    help="Fontsize of xticks. Default: matplotlib default")
# -------- End of Label and Tick options --------

# -------- Scaling options --------
parser.add_argument("--subplot_scale", action="store_true", dest="subplot_scale",
                    help="Scale feature x size by subplot x/y ratio. Default: off")
parser.add_argument("--track_group_scale", action="store_true", dest="track_group_scale",
                    help="Scale feature x size by track_group x/y ratio. Default: off")
# -------- End of Scaling options --------
# ---- End of Drawing options ----

# ---- Common options ----

# -------- Subplot adjustment options --------

parser.add_argument("--manual_figure_adjustment", action="store_true", dest="manual_figure_adjustment", default=False,
                    help="Adjust borders of figure manually using options below. Default: False, i.e. scaling is done automatically.")
parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float,
                    help="Adjust left border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_top", action="store", dest="subplots_adjust_top", type=float,
                    help="Adjust top border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_right", action="store", dest="subplots_adjust_right", type=float,
                    help="Adjust right border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_bottom", action="store", dest="subplots_adjust_bottom", type=float,
                    help="Adjust bottom border of subplots on the figure. Default: matplotlib defaults")
# -------- End of Subplot adjustment options --------

# -------- Figure size options --------
parser.add_argument("--figure_header_height", action="store", dest="figure_header_height",
                    type=float, default=0.0,
                    help="Height of figure header. Default: 0.0")
parser.add_argument("--figure_height_per_scaffold", action="store", dest="figure_height_per_scaffold",
                    type=float, default=0.5,
                    help="Height of figure per chromosome track. Default: 0.5")
parser.add_argument("--figure_width", action="store", dest="figure_width", type=float, default=10,
                    help="Width of figure in inches. Default: 10")
# -------- End of Figure size options --------

# -------- Miscellaneous options --------
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Print additional info to stdout")
# -------- End of Miscellaneous options --------

# ---- End of Common options ----
args = parser.parse_args()


try:
    feature_start_column_id = "start"
    feature_end_column_id = "end"
    feature_color_column_id = None

    if args.input_type in ["gtf", "gff"]:  # GTF and GFF formats. Only coordinates are used. Comment lines are allowed and must start from '#'
        feature_df = CollectionGFF(in_file=args.input, parsing_mode="only_coordinates")

    elif args.input_type in ["bed", "bedgraph"]:  # Bed format without track lines. All columns except first three are ignored. Comment lines are allowed and must start from '#'
        feature_df = CollectionBED(in_file=args.input, parsing_mode="coordinates_only", format="bed")

    elif args.input_type in ["bed_colored"]:  # Four column bed with following fields: "scaffold", start", "end", "color". Comment lines are not allowed.
        feature_df = CollectionBED(in_file=args.input, parsing_mode="complete", format="bed_colored")
        #feature_df.records.columns = pd.Index(["start", "end", "color"])
        #feature_df.records.index.name = "scaffold"
        feature_color_column_id = "color" if not args.enforce_default_color else None

    elif args.input_type in ["bed_track", "bedgraph_colored"]:  # Five-column bed with following fields: "scaffold", start", "end", "value", "color". The 'value' field is ignored. Comment lines are  NOT ALLOWED...
        feature_df = CollectionBED(in_file=args.input, parsing_mode="all", format="bed_track",)
        feature_color_column_id = "color" if not args.enforce_default_color else None

    elif args.input_type == "bed_with_header":  # 'BED' format without track lines, but WITH a header line. Not a true BED!
        feature_df = CollectionBED(in_file=args.input, parsing_mode="coordinates_only", format="bed",
                                   header_in_file=True)

    elif args.input_type == "table":  # A tab separated table WITH a header line. Scaffold, start and end columns must be set via script options. Comment lines are allowed and must start from '#'
        feature_df = CollectionBED(in_file=args.input, parsing_mode="coordinates_only", format="table",
                                   header_in_file=True,
                                   scaffold_column_name=args.scaffold_column_name,
                                   start_column_name=args.start_column_name,
                                   end_column_name=args.end_column_name,
                                   one_based_coordinates=args.one_based_coordinates)

    elif args.input_type == "table_colored":  # A tab separated table with a header line. Scaffold, start, end and color columns must be set via script options. Comment lines are  NOT ALLOWED.
        feature_df = CollectionBED(in_file=args.input, parsing_mode="colored", format="table",
                                   header_in_file=True,
                                   scaffold_column_name=args.scaffold_column_name,
                                   start_column_name=args.start_column_name,
                                   end_column_name=args.end_column_name,
                                   color_column_name=args.color_column_name,
                                   one_based_coordinates=args.one_based_coordinates
                                   )
        feature_color_column_id = "color" if not args.enforce_default_color else None

    elif args.input_type == "blast6":  # BLAST output 6 format
        feature_df = CollectionBLAST(in_file=args.input, parsing_mode="complete")
        feature_df.records.reset_index(level="query_id", inplace=True)
        feature_start_column_id = "target_start"
        feature_end_column_id = "target_end"

    elif args.input_type == "blast6_colored":  # BLAST output 6 format with header and additional color column
        feature_df = CollectionBLAST(in_file=args.input, parsing_mode="complete", format="tab6_colored", header=True)
        feature_df.records.reset_index(level="query_id", inplace=True)
        feature_start_column_id = "target_start"
        feature_end_column_id = "target_end"
        feature_color_column_id = args.color_column_name

    elif args.input_type == "STR": #
        feature_df = CollectionSTR(in_file=args.input, records=None, format="filtered_str", parsing_mode="all",
                                   black_list=(), white_list=())

        feature_df.records.set_index("scaffold_id", inplace=True)
        feature_df.records.index.name = "scaffold"

    else:
        raise ValueError(f"ERROR!!! Unrecognized input type ({args.input_type}). ")

except pd.errors.EmptyDataError:
    print("Empty input file. Silent exit.")   # try-except added to handle case when input file is empty without raising exception. For use in snakemake
    exit(0)

auxiliary_dict = Parsing.read_mace_auxiliary_input(len_file=args.scaffold_length_file,
                                                   whitelist_file=args.scaffold_whitelist,
                                                   max_scaffolds=args.max_scaffolds,
                                                   orderlist_file=args.scaffold_orderlist,
                                                   syn_file=args.scaffold_syn_file,
                                                   syn_file_key_column=args.syn_file_key_column,
                                                   syn_file_value_column=args.syn_file_value_column,
                                                   centromere_bed=args.centromere_bed,
                                                   highlight_bed=args.highlight_file,
                                                   legend_file=args.legend,
                                                   vert_track_group_file=None,
                                                   hor_track_group_file=None,
                                                   hor_track_subgroup_file=None)

records_df = Parsing.resolve_mace_single_genome_input(auxiliary_dict, records_df=feature_df.records)

Visualization.draw_features({"features": records_df}, auxiliary_dict["len_df"],
                            auxiliary_dict["orderlist_series"],
                            args.output_prefix,
                            legend=Visualization.feature_legend(auxiliary_dict["legend_df"], colormap=args.legend_colormap),
                            #legend_df=legend_df,
                            centromere_df=auxiliary_dict["centromere_df"],
                            highlight_df=auxiliary_dict["highlight_df"],
                            figure_width=args.figure_width,
                            figure_height_per_scaffold=0.5,
                            figure_header_height=args.figure_header_height,
                            dpi=300,
                            #colormap=None, thresholds=None, colors=None, background=None,
                            default_color=args.default_color,
                            title=args.title,
                            extensions=args.output_formats,
                            feature_shape=args.feature_shape,
                            feature_start_column_id=feature_start_column_id,
                            feature_end_column_id=feature_end_column_id,
                            feature_color_column_id=feature_color_column_id,
                            feature_length_column_id="length",
                            subplots_adjust_left=args.subplots_adjust_left,
                            subplots_adjust_bottom=args.subplots_adjust_bottom,
                            subplots_adjust_right=args.subplots_adjust_right,
                            subplots_adjust_top=args.subplots_adjust_top,
                            show_track_label=not args.hide_track_label,
                            show_trackgroup_label=True,
                            subplot_scale=args.subplot_scale,
                            track_group_scale=args.track_group_scale,
                            stranded_tracks=args.stranded,
                            rounded_tracks=args.rounded,
                            stranded_end_tracks=args.stranded_end,
                            fill_empty_tracks=args.fill_empty_tracks,
                            empty_color=args.empty_color,
                            xtick_fontsize=args.x_tick_fontsize,
                            subplot_title_fontsize=args.title_fontsize,
                            subplot_title_fontweight='bold',
                            x_tick_type=args.x_tick_type,
                            autoscale_figure=False if args.manual_figure_adjustment else True,
                            )
