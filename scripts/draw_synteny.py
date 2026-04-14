#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from copy import deepcopy
from pathlib import Path
from collections import OrderedDict

import pandas as pd

from distinctipy import distinctipy

from RouToolPa.Parsers.PSL import CollectionPSL
from RouToolPa.Parsers.BED import CollectionBED

from MACE.Routines import Parsing
from MACE.Routines import Synteny
from MACE.Routines import Visualization


def bed_dict_to_xlsx(bed_dict, genome_auxiliary_dict, output_prefix):
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
        for scaffold in genome_auxiliary_dict[species]["scaffold_color_df"].index:
            species_format_dict[species][scaffold] = workbook.add_format(
                {'bg_color': genome_auxiliary_dict[species]["scaffold_color_df"].loc[scaffold, "color"]})

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
        #workbook.get_worksheet_by_name('species').autofit()  # NOT verified
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


parser = argparse.ArgumentParser()
# ---- Input/Output options ----

# -------- Main Input options --------
parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True,
                    help="Input directory with data. Must contain subfolder for each genome. "
                         "Subfolders should have same name as genomes in --genome_orderlist"
                         "Each subfolder should contain: *.whitelist, *.len and synteny file "
                         "(except for the reference genome). *.orderlist, *.invertlist and *.syn file are optional."
                         "The folder of the reference genome (set via --reference) may additionally include centromere.bed")
parser.add_argument("--synteny_format", action="store", dest="synteny_format", default="psl",
                    help="Format of the synteny file. Allowed: psl(default), bed, bed_with_color")

# -------- Output options -------
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=("png", "svg"),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg,png")
# -------- End of Output options ------

# ---- End of Input/Output options ----

# -------- Reference-Related options --------
parser.add_argument("--reference", action="store", dest="reference", required=True, type=str,
                    help="Name of the  reference genome")
#parser.add_argument("--reference_centromere_bed", action="store", dest="reference_centromere_bed", required=False,
#                    type=str,
#                    help="Bed file with coordinates of centromeres in reference")

# -------- End of Reference-Related options --------

# -------- Query-Related options --------
parser.add_argument("--query_orderlist", action="store", dest="query_orderlist", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of labels for query genomes.")
parser.add_argument("--invert_genome_order", action="store_true", dest="invert_genome_order", default=False,
                    help="Invert order of the genomes in the --genome_orderlist. Default: False")

parser.add_argument("--use_original_colors", action="store_true", dest="use_original_colors", default=False,
                    help="Use colors from .color file. If not set colors will be select ")
parser.add_argument("--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in synonym file. Default: 0")
parser.add_argument("--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in synonym file synonym. Default: 1")
# --------End of Query-Related options --------

# ---- Data Processing options --------
# -------- Synteny Block Filtering and Processing options --------

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
# -------- End of Synteny Block Filtering and Processing options --------
parser.add_argument("--invert_coordinates_for_target_negative_strand", action="store_true",
                    dest="invert_coordinates_for_target_negative_strand",
                    default=False,
                    help="Invert coordinates for target negative strand. "
                         "Set only if input PSL file doesn't follow specification (target should have an inverted coordinate "
                         "system for negative synteny blocks). Default: False")

# ---- End of Data Processing options --------
# ---- Drawing options ----

# -------- Chromosome track options --------
parser.add_argument("--stranded", action="store_true", dest="stranded", default=False,
                    help="Stranded features and tracks. Default: False")
parser.add_argument("--rounded", action="store_true", dest="rounded", default=False,
                    help="Rounded tracks. Default: False")
parser.add_argument("--stranded_end", action="store_true", dest="stranded_end", default=False,
                    help="Stranded ends for tracks. Works only if --stranded is set. Default: False")
parser.add_argument("--chromosome_height", action="store", dest="chromosome_height", default=9, type=float,
                    help="Height of chromosomes on the plot. Increase or decrease this parameter to make chromosomes "
                         "thicker or thinner. Default: 9")
# -------- End of Chromosome track options --------

# -------- Legend options --------
parser.add_argument("--hide_legend", action="store_true", dest="hide_legend", default=False,
                    help="Don't draw legend. Default: False")
# -------- End of Legend options --------

# -------- Title options --------
parser.add_argument("-l", "--title", action="store", dest="title", default="Synteny",
                    help="Suptitle of figure. Default: 'Synteny'")
parser.add_argument("--title_fontsize", action="store", dest="title_fontsize", default=20, type=int,
                    help="Fontsize of the figure. Default: 20")
# -------- End of Title options --------

# -------- Label and Tick options --------
parser.add_argument("--hide_track_label", action="store_true", dest="hide_track_label", default=False,
                    help="Hide track label. Default: False")
parser.add_argument("--x_tick_fontsize", action="store", dest="x_tick_fontsize", type=int, default=None,
                    help="Fontsize of xticks. Default: matplotlib default")
# -------- End of Label and Tick options --------

# -------- Scaling options --------
parser.add_argument("--subplot_scale", action="store_true", dest="subplot_scale",
                    help="Scale feature x size by subplot x/y ratio. Default: off")
parser.add_argument("--track_group_scale", action="store_true", dest="track_group_scale",
                    help="Scale feature x size by track_group x/y ratio. Default: off")
parser.add_argument("--ymax_multiplier", action="store", dest="ymax_multiplier", type=float, default=1.0,
                    help="Multiplier for y max limit of figure. Default: 1.0")
# -------- End of Scaling options --------

# ---- End of Drawing options --

# ---- Common options ----
# -------- Subplot adjustment options --------
parser.add_argument("--manual_figure_adjustment", action="store_true", dest="manual_figure_adjustment", default=False,
                    help="Adjust borders of figure manually using options below. Default: False, i.e. scaling is done automatically.")
parser.add_argument("--subplots_adjust_left", action="store", dest="subplots_adjust_left", type=float, default=0.2,
                    help="Adjust left border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_top", action="store", dest="subplots_adjust_top", type=float,
                    help="Adjust top border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_right", action="store", dest="subplots_adjust_right", type=float,
                    help="Adjust right border of subplots on the figure. Default: matplotlib defaults")
parser.add_argument("--subplots_adjust_bottom", action="store", dest="subplots_adjust_bottom", type=float,
                    help="Adjust bottom border of subplots on the figure. Default: matplotlib defaults")
# -------- End of Subplot adjustment options --------

# -------- Figure size options --------

# -------- End of Figure size options --------
parser.add_argument("--figure_header_height", action="store", dest="figure_header_height", type=int, default=0,
                    help="Additional height of figure to account for header. Default: 0")
parser.add_argument("--figure_height_per_scaffold", action="store", dest="figure_height_per_scaffold",
                    type=float, default=0.5,
                    help="Height of figure per chromosome track. Default: 0.5")
parser.add_argument("--figure_width", action="store", dest="figure_width", type=float, default=15,
                    help="Width of figure in inches. Default: 15")
# ---- End of Common options ----

# -------- Miscellaneous options --------
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Print additional info to stdout")
# -------- End of Miscellaneous options --------

args = parser.parse_args()

reference = args.reference
query_list = args.query_orderlist

data_dir = args.input_dir
data_dir_path = Path(data_dir)

genome_list = query_list + [reference]
syn_file_key_column, syn_file_value_column = args.syn_file_key_column, args.syn_file_value_column

synteny_format = args.synteny_format
if args.verbose:
    for genome in genome_list:
        print(data_dir_path / genome)

genome_auxiliary_dict = {genome: Parsing.read_mace_auxiliary_input(
                                              len_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".len", ".fasta", ".fai"]),
                                              whitelist_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".whitelist"]),
                                              orderlist_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".orderlist"]),
                                              invertlist_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".invertlist"]),
                                              syn_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".syn"]),
                                              syn_file_key_column=args.syn_file_key_column,
                                              syn_file_value_column=args.syn_file_value_column,
                                              centromere_bed=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".centromere.bed"]),
                                              highlight_bed=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".highlight.bed"]),
                                              scaffold_color_file=Parsing.get_filenames_for_extension(data_dir_path / genome, extension_list=[".scaffold.color"]),
                                              vert_track_group_file=None,
                                              hor_track_group_file=None,
                                              hor_track_subgroup_file=None) for genome in genome_list}

if not args.use_original_colors:
    color_number = max([len(genome_auxiliary_dict[query]["orderlist_series"]) for query in query_list])
    colors = distinctipy.get_colors(color_number)
    color_list = list(map(Visualization.rgb_tuple_to_hex, colors))

    for species in query_list:
        genome_auxiliary_dict[species]["scaffold_color_df"] = pd.DataFrame()

        genome_auxiliary_dict[species]["scaffold_color_df"]["scaffold"] = genome_auxiliary_dict[species]["orderlist_series"]
        genome_auxiliary_dict[species]["scaffold_color_df"]["color"] = color_list[:len(genome_auxiliary_dict[species]["orderlist_series"])]
        genome_auxiliary_dict[species]["scaffold_color_df"].set_index("scaffold", inplace=True)
        genome_auxiliary_dict[species]["scaffold_color_df"].to_csv("{}.{}.chr_colors.tsv".format(args.output_prefix, species), sep="\t", header=True, index=True)

bed_col_dict = OrderedDict()

if args.synteny_format == "psl":
    for query in query_list:
        print(Parsing.get_filenames_for_extension(data_dir_path / query, extension_list=["psl", "psl.gz"]),)
    psl_col_dict = {query: CollectionPSL(in_file=Parsing.get_filenames_for_extension(data_dir_path / query,
                                                                                     extension_list=["psl", "psl.gz"]),
                                         parsing_mode="coordinates_only",
                                         #target_syn_dict=syn_df_dict[reference].to_dict(),
                                         target_black_list=None,
                                         target_white_list=genome_auxiliary_dict[reference]["whitelist_series"],
                                         #query_syn_dict=syn_df_dict[query].to_dict(),
                                         query_black_list=None,
                                         query_white_list=genome_auxiliary_dict[query]["whitelist_series"],
                                         invert_coordinates_for_target_negative_strand=args.invert_coordinates_for_target_negative_strand
                                         ) for query in query_list}
    for query in query_list:
        psl_col_dict[query].records["tName"] = psl_col_dict[query].records["tName"].replace(genome_auxiliary_dict[reference]["syn_dict"])
        psl_col_dict[query].records["qName"] = psl_col_dict[query].records["qName"].replace(genome_auxiliary_dict[query]["syn_dict"])

    for species in query_list:
        bed_col_dict[species] = CollectionBED(
            records=psl_col_dict[species].records[["tName", "tStart", "tEnd", "qName", "qStart", "qEnd", "strand"]],
            records_columns=["scaffold", "start", "end", "query", "query_start", "query_end", "strand"],
            format="bed_synteny_track",
            parsing_mode="all")
        bed_col_dict[species].records.set_index("scaffold", inplace=True)

        bed_col_dict[species].records["color"] = bed_col_dict[species].records["query"].replace(genome_auxiliary_dict[species]["scaffold_color_df"]["color"])
        bed_col_dict[species].records = bed_col_dict[species].records.sort_values(by=["scaffold", "start", "end", ])

elif args.synteny_format in ["bed", "bed_with_color"]: # TODO: REWRITE!!!!!!!!!!!!!!
    bed_file_dict = {query: Parsing.expand_path(bed, skip=not args.expand_paths) for query, bed in zip(query_list, args.input)}

    for species in query_list:
        bed_col_dict[species] = CollectionBED(in_file=Parsing.get_filenames_for_extension(data_dir_path / species,
                                                                                          extension_list=["bed", "bed.gz"]), header_in_file=True,
                                              format="bed_synteny_track", parsing_mode="all",
                                              scaffold_syn_dict=genome_auxiliary_dict[reference]["syn_dict"],
                                              rename_dict={"query": genome_auxiliary_dict[species]["syn_dict"]} if genome_auxiliary_dict[species]["syn_dict"] is not None else None)

        bed_col_dict[species].records = bed_col_dict[species].records.sort_values(by=["scaffold", "start", "end", ])
        if args.synteny_format != "bed_with_color":
            bed_col_dict[species].records["color"] = bed_col_dict[species].records["query"].replace(genome_auxiliary_dict[species]["scaffold_color_df"]["color"])
else:
    raise ValueError("ERROR!!! Unrecognized format of the input file(s)!")

for genome in genome_list:
    Parsing.resolve_mace_single_genome_input(genome_auxiliary_dict[genome])

query_scaffold_id_column_name = "query"
query_start_column_name = "query_start"
query_end_column_name = "query_end"
strand_column_name = "strand"

for species in query_list:
    print("Inverting (if necessary) {0} scaffolds...".format(species))
    print("Inverting query coordinates in synteny file...")

    bed_col_dict[species].records = Synteny.invert_coordinates_in_synteny_table(bed_col_dict[species].records,
                                                                                genome_auxiliary_dict[species]["invertlist_series"],
                                                                                genome_auxiliary_dict[species]["preinvert_len_df"],
                                                                                query_scaffold_id_column_name,
                                                                                query_start_column_name,
                                                                                query_end_column_name,
                                                                                strand_column_name,
                                                                                genome_auxiliary_dict[species]["inverted_scaffold_label"])
#--------------------------------------------------------------------------------------

query_species_color_df_dict = {sp: genome_auxiliary_dict[sp]["scaffold_color_df"] for sp in query_list}

common_visualization_options_dict = {"legend": None if args.hide_legend else Visualization.chromosome_legend(query_species_color_df_dict,
                                                                                                             genome_auxiliary_dict[reference]["orderlist_series"]),
                                     "centromere_df": genome_auxiliary_dict[reference]["centromere_df"],
                                     "highlight_df": genome_auxiliary_dict[reference]["highlight_df"],
                                     "figure_width": args.figure_width,
                                     "figure_height_per_scaffold": args.figure_height_per_scaffold,
                                     "dpi": 300,
                                     "default_color": "red",  # TODO: check if it is possible to remove it
                                     "title": args.title,
                                     "extensions": args.output_formats,
                                     "feature_start_column_id": "start",
                                     "feature_end_column_id": "end",
                                     "feature_color_column_id": "color",
                                     "feature_length_column_id": "length",
                                     "feature_height_fraction": 0.7,
                                     "subplots_adjust_left": args.subplots_adjust_left,
                                     "subplots_adjust_bottom": args.subplots_adjust_bottom,
                                     "subplots_adjust_right": args.subplots_adjust_right,
                                     "subplots_adjust_top": args.subplots_adjust_top,
                                     "autoscale_figure": False if args.manual_figure_adjustment else True,
                                     "show_track_label": not args.hide_track_label,
                                     "show_trackgroup_label": True,
                                     "close_figure": True,
                                     "subplot_scale": False,
                                     "track_group_scale": False,
                                     "track_group_distance": 2,
                                     "xmax_multiplier": 1.3,
                                     "ymax_multiplier": args.ymax_multiplier,
                                     "figure_header_height": args.figure_header_height,
                                     "rounded_tracks": args.rounded,
                                     "xtick_fontsize": args.x_tick_fontsize,
                                     "subplot_title_fontsize": args.title_fontsize,
                                     "subplot_title_fontweight": 'bold',}

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
                                genome_auxiliary_dict[reference]["len_df"],  #reference_scaffold_length_df,
                                genome_auxiliary_dict[reference]["orderlist_series"],
                                "{0}.initial_min_block_len_{1}".format(args.output_prefix, min_block_length),
                                stranded_tracks=args.stranded,
                                stranded_end_tracks=args.stranded_end,
                                **common_visualization_options_dict
                                )
    bed_dict_to_xlsx(prefiltered_bed_col_dict, genome_auxiliary_dict,
                     '{0}.initial_min_block_len_{1}'.format(args.output_prefix, min_block_length))

    for secondary_min_block_len in args.secondary_min_block_len_list:
        for max_dist_between_short_blocks in args.max_dist_between_short_blocks_list:
            filtered_bed_col_dict = deepcopy(prefiltered_bed_col_dict)

            for species in filtered_bed_col_dict:
                filtered_bed_col_dict[species].records = Synteny.filter_isolated_short_blocks(filtered_bed_col_dict[species].records,
                                                                                              min_block_len=secondary_min_block_len,
                                                                                              max_dist_between_short_blocks=max_dist_between_short_blocks)

            second_stage_output_suffix = "initial_min_block_len_{0}.secondary_min_block_len_{1}.max_dist_between_short_blocks_{2}".format(min_block_length,
                                                                                                                                          secondary_min_block_len,
                                                                                                                                          max_dist_between_short_blocks)

            for species in filtered_bed_col_dict:  # bed_col_dict:
                # ---- save original blocks to bed ----
                filtered_bed_col_dict[species].records.to_csv("{0}.{1}.to.{2}.{3}.tsv".format(args.output_prefix,
                                                                                              species,
                                                                                              reference,
                                                                                              second_stage_output_suffix),
                                                              sep="\t",
                                                              header=True, index=True)
                # -----

            Visualization.draw_features(filtered_bed_col_dict,
                                        genome_auxiliary_dict[reference]["len_df"],
                                        genome_auxiliary_dict[reference]["orderlist_series"],
                                        "{0}.{1}".format(args.output_prefix, second_stage_output_suffix),
                                        stranded_tracks=args.stranded,
                                        stranded_end_tracks=args.stranded_end,
                                        **common_visualization_options_dict
                                        )
            bed_dict_to_xlsx(filtered_bed_col_dict, genome_auxiliary_dict,
                             "{0}.{1}".format(args.output_prefix,
                                              second_stage_output_suffix))

            for max_dist_between_blocks in args.max_dist_between_blocks_list:
                for species in filtered_bed_col_dict:
                    filtered_bed_col_dict[species].records = Synteny.merge_adjacent_blocks(filtered_bed_col_dict[species].records,
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
                                            genome_auxiliary_dict[reference]["len_df"],
                                            genome_auxiliary_dict[reference]["orderlist_series"],
                                            "{0}.{1}.{2}".format(args.output_prefix,
                                                                 second_stage_output_suffix,
                                                                 third_stage_output_suffix),
                                            stranded_tracks=False,
                                            stranded_end_tracks=False,
                                            **common_visualization_options_dict
                                            )
                bed_dict_to_xlsx(filtered_bed_col_dict, genome_auxiliary_dict,
                                 "{0}.{1}.{2}".format(args.output_prefix, second_stage_output_suffix, third_stage_output_suffix))

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
                                                genome_auxiliary_dict[reference]["len_df"],
                                                genome_auxiliary_dict[reference]["orderlist_series"],
                                                "{0}.{1}.{2}.{3}".format(args.output_prefix,
                                                                         second_stage_output_suffix,
                                                                         third_stage_output_suffix,
                                                                         forth_stage_output_suffix),
                                                stranded_tracks=False,
                                                stranded_end_tracks=False,
                                                **common_visualization_options_dict
                                                )
                    bed_dict_to_xlsx(filtered_bed_col_dict, genome_auxiliary_dict,
                                     "{0}.{1}.{2}.{3}".format(args.output_prefix,
                                                              second_stage_output_suffix,
                                                              third_stage_output_suffix,
                                                              forth_stage_output_suffix))
