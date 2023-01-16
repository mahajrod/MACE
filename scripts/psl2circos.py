#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os

import shutil
import argparse
from collections import OrderedDict
import importlib.resources as resources

from pathlib import Path
import pandas as pd
import numpy as np

from distinctipy import distinctipy
from functools import partial
from RouToolPa.Collections.General import SynDict, IdList
from MACE.Routines import Circos

from RouToolPa.Parsers.PSL import CollectionPSL
from RouToolPa.Parsers.BED import CollectionBED


def parse_id_list(string):
    return pd.read_csv(string, header=None, squeeze=True) if Path(string).exists() else pd.Series(string.split(","))


def split_comma_separated_list(string):
    return string.split(",")


def rgb_tuple_to_hex(rgb_tuple):
    color_code = "#"
    for i in [0, 1, 2]:
        color_code += "{:02X}".format(int(255 * rgb_tuple[i]))

    return color_code


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with synteny data")
parser.add_argument("-f", "--format", action="store", dest="format", default="psl",
                    help="Format of input file. Allowed: psl(default), bed.")

parser.add_argument("--query_label", action="store", dest="query_label", required=True,
                    help="Labels for query genomes.")
parser.add_argument("--reference_label", action="store", dest="reference_label", required=True, type=str,
                    help="Label of reference genome")

parser.add_argument("--query_white_list", action="store", dest="query_white_list", default=None,
                    type=parse_id_list,
                    help="Comma-separated list or file containing the only scaffolds from query genomes to include."
                         "Default: all")
parser.add_argument("--query_black_list", action="store", dest="query_black_list", default=None,
                    type=parse_id_list,
                    help="Comma-separated list or file containing scaffolds from query genome to exclude. "
                         "Default: not set")

parser.add_argument("--target_white_list", action="store", dest="target_white_list", default=None,
                    type=parse_id_list,
                    help="Comma-separated list or file containing scaffolds from target genome to include. "
                         "Default: all")
parser.add_argument("--target_black_list", action="store", dest="target_black_list",
                    default=None,
                    type=parse_id_list,
                    help="Comma-separated list or file containing scaffolds from target genome to exclude. "
                         "Default: not set")

parser.add_argument("--query_syn_file", action="store", dest="query_syn_file",
                    help="Files with scaffold id synonyms for query genomes")
parser.add_argument("--target_syn_file", action="store", dest="target_syn_file",
                    help="File with scaffold id synonyms for target genome")
parser.add_argument("--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in synonym file. Default: 0")
parser.add_argument("--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in synonym file synonym. Default: 1")

parser.add_argument("-z", "--target_scaffold_order_list", action="store", dest="target_scaffold_order_list",
                    required=True,
                    type=parse_id_list,
                    help="Comma-separated list of scaffolds from target genome to draw first and exactly in same order."
                         "Default: not set")
parser.add_argument("--query_scaffold_order_list", action="store", dest="query_scaffold_order_lists", required=True,
                    type=parse_id_list,
                    help="Comma-separated list of scaffolds from query genome to draw first and exactly in same order."
                         "Default: not set")
parser.add_argument("-n", "--reference_scaffold_length_file", action="store", dest="reference_scaffold_length_file",
                    required=True,
                    help="File with lengths of scaffolds for reference genome")
parser.add_argument("--query_scaffold_length_file", action="store", dest="query_scaffold_length_file",
                    required=True,
                    help="Comma-separated list of labels for query genome. Must follow the same order as for bed files")

parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    type=lambda s: Path(s),
                    help="Path to output directory")

parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=split_comma_separated_list,
                    default=("png", ),
                    help="Comma-separated list of formats (supported by matlotlib) of "
                         "output figure.Default: svg,png")

parser.add_argument("-l", "--title", action="store", dest="title", default="Coverage",
                    help="Suptitle of figure. Default: 'Coverage'")


args = parser.parse_args()

if not args.output_dir.exists():
    args.output_dir.mkdir()

target_syn_dict = SynDict(filename=args.target_syn_file,
                          key_index=args.syn_file_key_column,
                          value_index=args.syn_file_value_column) if args.target_syn_file is not None else {}
query_syn_dict = SynDict(filename=args.query_syn_file,
                         key_index=args.syn_file_key_column,
                         value_index=args.syn_file_value_column) if args.query_syn_file is not None else {}

if args.format == "psl":
    collection = CollectionPSL(in_file=args.input, parsing_mode="coordinates_only",
                               target_syn_dict=target_syn_dict,
                               target_black_list=args.target_black_list,
                               target_white_list=args.target_white_list,
                               query_syn_dict=query_syn_dict,
                               query_black_list=args.query_black_list,
                               query_white_list=args.query_white_list
                               )
elif args.format == "bed":
    collection = CollectionBED(in_file=args.input, header_in_file=True,
                               format="bed_synteny_track", parsing_mode="all",
                               scaffold_white_list=args.target_white_list,
                               scaffold_black_list=args.target_black_list,
                               scaffold_syn_dict=target_syn_dict,
                               rename_dict={"query": query_syn_dict},
                               white_list_dict={"query": args.query_white_list},
                               black_list_dict={"query": args.query_black_list}
                               )
else:
    raise ValueError("ERROR!!! Unknown collection type")

links_df = Circos.convert_bed_synteny_track_to_links(collection, type=args.format)
links_df.to_csv(args.output_dir / "links")



"""
reference = args.reference_label
query_list = args.query_labels
interquery_list = args.interquery_bed_labels
species_chr_syn_dict = {query: pd.read_csv(syn_file,
                                           usecols=(args.syn_file_key_column,args.syn_file_value_column),
                                           header=None,
                                           sep="\t", squeeze=True, index_col=args.syn_file_key_column) if syn_file else None for query, syn_file in zip(query_list, args.query_scaffold_syn_files)}

species_chr_syn_dict[reference] = pd.read_csv(args.reference_scaffold_syn_file,
                                              usecols=(args.syn_file_key_column,args.syn_file_value_column,),
                                              header=None,
                                              sep="\t", squeeze=True, index_col=args.syn_file_key_column) if args.reference_scaffold_syn_file else None

if args.query_scaffold_white_lists:
    species_white_list_dict = {query: pd.read_csv(white_list_file, header=None, squeeze=True) if white_list_file else None for query, white_list_file in zip(query_list, args.query_scaffold_white_lists)}
else:
    species_white_list_dict = {query: None for query in query_list}

species_white_list_dict[reference] = args.reference_scaffold_white_list #pd.read_csv(args.reference_scaffold_white_list, header=None, squeeze=True) if args.reference_scaffold_white_list else None

if args.query_scaffold_black_lists:
    species_black_list_dict = {query: pd.read_csv(black_list_file, header=None, squeeze=True) if black_list_file else None for query, black_list_file in zip(query_list, args.query_scaffold_black_lists)}
else:
    species_black_list_dict = {query: None for query in query_list}

species_black_list_dict[reference] = args.reference_scaffold_black_list #pd.read_csv(args.reference_scaffold_black_list, header=None, squeeze=True) if args.reference_scaffold_black_list else None

print(args.input)
for query in query_list + [reference]:
    print(species_chr_syn_dict[query])

#psl_col_dict = {query: CollectionPSL(in_file=psl, parsing_mode="coordinates_only",
#                                     target_syn_dict=species_chr_syn_dict[reference],
#                                     target_black_list=species_black_list_dict[reference],
#                                     target_white_list=species_white_list_dict[reference],
#                                     query_syn_dict=species_chr_syn_dict[query],
#                                     query_black_list=species_black_list_dict[query],
#                                     query_white_list=species_white_list_dict[query]
#                                     ) for query, psl in zip(query_list, args.input)}

species_orderlist_dict = {query: pd.read_csv(orderlist_file, header=None, squeeze=True).iloc[::-1] for query, orderlist_file in zip(query_list, args.query_scaffold_order_lists)}
species_orderlist_dict[reference] = args.reference_scaffold_order_list[::-1]

color_number = max([len(species_orderlist_dict[query]) for query in query_list])
colors = distinctipy.get_colors(color_number)
color_list = list(map(rgb_tuple_to_hex, colors))


species_color_df_dict = OrderedDict()
for species in query_list:
    species_color_df_dict[species] = pd.DataFrame()
    species_color_df_dict[species]["scaffold"] = species_orderlist_dict[species]
    species_color_df_dict[species]["color"] = color_list[:len(species_orderlist_dict[species])]
    species_color_df_dict[species].set_index("scaffold", inplace=True)

# species_color_df_dict[species]

reference_scaffold_length_df = pd.read_csv(args.reference_scaffold_length_file, sep='\t',
                                           header=None, names=("scaffold", "length"), index_col=0)
reference_scaffold_length_df.index = pd.Index(list(map(str, reference_scaffold_length_df.index)))
reference_scaffold_length_df.rename(index=species_chr_syn_dict[reference], inplace=True)

scaffold_length_df_dict = OrderedDict()
scaffold_length_df_dict[reference] = reference_scaffold_length_df
for query, length_file in zip(query_list, args.query_scaffold_length_files):
    scaffold_length_df_dict[query] = pd.read_csv(args.query_scaffold_length_files, sep='\t',
                                           header=None, names=("scaffold", "length"), index_col=0)
    scaffold_length_df_dict[query].index = pd.Index(list(map(str, reference_scaffold_length_df.index)))
    scaffold_length_df_dict[query].rename(index=species_chr_syn_dict[query], inplace=True)

bed_col_dict = OrderedDict()

bed_file_dict = {query: bed for query, bed in zip(query_list, args.input)}

for species in query_list:
    bed_col_dict[species] = CollectionBED(in_file=bed_file_dict[species], header_in_file=True,
                                          format="bed_synteny_track", parsing_mode="all",
                                          scaffold_syn_dict=species_chr_syn_dict[reference] if species_chr_syn_dict[reference] is not None else None,
                                          rename_dict={"query": species_chr_syn_dict[species]} if species_chr_syn_dict[species] is not None else None)
    #bed_col_dict[species].records.set_index("scaffold", inplace=True)

    #bed_col_dict[species].records["color"] = bed_col_dict[species].records["query"].replace(species_color_df_dict[species]["color"])
    print(bed_col_dict[species].records)
    bed_col_dict[species].records.sort_values(by=["scaffold", "start", "end", ], inplace=True)
    # bed_col_dict[species].records["color"][~bed_col_dict[species].records["color"].in(color_list)] = ""

bed_col_dict["interquery"] = CollectionBED(in_file=bed_file_dict[args.interquery_bed], header_in_file=True,
                                          format="bed_synteny_track", parsing_mode="all",
                                          scaffold_syn_dict=species_chr_syn_dict[interquery_list[0]] if species_chr_syn_dict[interquery_list[0]] is not None else None,
                                          rename_dict={"query": species_chr_syn_dict[interquery_list[1]]} if species_chr_syn_dict[interquery_list[1]] is not None else None)
print(bed_col_dict["interquery"].records)
bed_col_dict["interquery"].records.sort_values(by=["scaffold", "start", "end", ], inplace=True)

query_species_color_df_dict = {sp: species_color_df_dict[sp] for sp in query_list}

for query in bed_col_dict:
    print(bed_col_dict[query].records)


links_df_dict = {species: Circos.convert_bed_synteny_track_to_links(bed_col_dict[species]) for species in bed_col_dict}

karyotype_df = pd.concat([Circos.convert_length_df_to_karyotype(scaffold_length_df_dict[species]) for species in scaffold_length_df_dict],
                     ignore_index=True)

output_dir = Path("circos").absolute()

with resources.path("MACE.resources.circos.syntheny1", "bands.conf") as bands, \
     resources.path("MACE.resources.circos.syntheny1", "ideogram.conf") as ideogram, \
     resources.path("MACE.resources.circos.syntheny1", "ideogram.label.conf") as ideogram_label, \
     resources.path("MACE.resources.circos.syntheny1", "ideogram.position.conf") as ideogram_position, \
     resources.path("MACE.resources.circos.syntheny1", "ticks.conf") as ticks:
    resource_dict = {"bands": bands,
                     "ideogram":  ideogram,
                     "ideogram_label": ideogram_label,
                     "ideogram_position": ideogram_position,
                     "ticks": ticks}

for scaffold in species_orderlist_dict[reference]:
    scaffold_dir = output_dir / scaffold

    links_path_dict = {species: scaffold_dir / "species.links" for species in links_df_dict}

    Circos.safe_mkdir(scaffold_dir)
    query_scaffolds_dict = OrderedDict()
    for species in query_list:
        links_df_dict[species][links_df_dict[species].index == reference].write(links_df_dict[species], sep="\t", header=False, index=False)
        query_scaffolds_dict[species] = links_df_dict[species].index[links_df_dict[species].index == reference]
    query_scaffolds_dict[]
    for resource in resource_dict:
        shutil.copy(resource_dict[resource], scaffold_dir)

    aligned_query_scaffolds = links_df_dict[links_df.index == scaffold]["query"]
    aligned_ref_query_scaffolds = aligned_query_scaffolds[aligned_query_scaffolds.isin(species_orderlist_dict[interquery_list[0]])]
    links_df[(links_df.index == scaffold) & (links_df.index.isin(aligned_ref_query_scaffolds))]

for species in bed_col_dict:
    bed_col_dict[species].records.to_csv("{0}.{1}.to.{2}.tsv".format(args.output_prefix, species, reference), sep="\t",
                                         header=True, index=True)


"""