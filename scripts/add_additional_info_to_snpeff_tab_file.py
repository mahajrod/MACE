#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from MACE.General.GeneralCollections import SynDict
from MACE.Parsers.SNPeff import CollectionSNPeff


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_snpeff", action="store", dest="input_snpeff", required=True,
                    help="Input snpeff tab file")
parser.add_argument("-o", "--output_snpeff", action="store", dest="output_snpeff", required=True,
                    help="Output snpeff file")

parser.add_argument("-a", "--alias_file", action="store", dest="alias_file",
                    help="File with gene name  aliases")
parser.add_argument("--alias_key", action="store", dest="alias_key", type=int, default=0,
                    help="Column with gene names in alias file (0-based). Default - 0")
parser.add_argument("--alias_value", action="store", dest="alias_value", type=int, default=1,
                    help="Column with gene aliases in alias file (0-based). Default - 1")

parser.add_argument("-f", "--function_file", action="store", dest="function_file",
                    help="File with gene functions")
parser.add_argument("--function_key", action="store", dest="function_key", type=int, default=0,
                    help="Column with gene names in function file (0-based). Default - 0")
parser.add_argument("--function_value", action="store", dest="function_value", type=int, default=1,
                    help="Column with gene functions in function file (0-based). Default - 1")

parser.add_argument("-d", "--description_file", action="store", dest="description_file",
                    help="File with descriptions")
parser.add_argument("--description_key", action="store", dest="description_key", type=int, default=0,
                    help="Column with gene names in description file (0-based). Default - 0")
parser.add_argument("--description_value", action="store", dest="description_value", type=int, default=1,
                    help="Column with description in description file (0-based). Default - 1")

parser.add_argument("-b", "--biochemical_file", action="store", dest="biochemical_file",
                    help="File with biochemical pathways")
parser.add_argument("--biochemical_key", action="store", dest="biochemical_key", type=int, default=0,
                    help="Column with gene names in biochemical pathway file (0-based). Default - 0")
parser.add_argument("--biochemical_value", action="store", dest="biochemical_value", type=int, default=1,
                    help="Column with biochemical pathway in biochemical pathway file (0-based). Default - 1")

args = parser.parse_args()

snpeff_collection = CollectionSNPeff(input_file=args.input_snpeff, from_file=True, filetype='tab')

if args.alias_file:
    alias_dict = SynDict()
    alias_dict.read(args.alias_file, split_values=True, key_index=args.alias_key, value_index=args.alias_value)
    snpeff_collection.add_gene_name_aliases(alias_dict)

if args.function_file:
    function_dict = SynDict()
    function_dict.read(args.function_file, split_values=True, key_index=args.function_key,
                       value_index=args.function_value)
    snpeff_collection.add_gene_functions(function_dict)

if args.description_file:
    description_dict = SynDict()
    description_dict.read(args.description_file, split_values=True, key_index=args.description_key,
                          value_index=args.description_value)
    snpeff_collection.add_gene_description(description_dict)

if args.biochemical_file:
    biochemical_dict = SynDict()
    biochemical_dict.read(args.biochemical_file, split_values=True, key_index=args.biochemical_key,
                          value_index=args.biochemical_value, allow_repeats_of_key=True)
    snpeff_collection.add_biochemical_pathway(biochemical_dict)

snpeff_collection.write(args.output_snpeff)



