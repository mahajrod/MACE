#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from Bio import SeqIO
from BCBio import GFF

from MACE.Parsers.VCF import CollectionVCF, ReferenceGenome


def list_from_str(s):
    return s.split(",")


def figsize_from_str(s):
    try:
        x, y = map(int, s.split(','))
        return x, y
    except:
        raise argparse.ArgumentTypeError("Figsize must be x,y")


def parse_synonyms(s):
    if not s:
        return None
    temp = s.split(",")
    synonyms_dict = {}
    for entry in temp:
        t = entry.split(":")
        synonyms_dict[t[0]] = t[1]
    return synonyms_dict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with mutations.")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", default="variant_location_pie",
                    help="Prefix of output file with pie chart")
parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=300,
                    help="Dpi of figure")
parser.add_argument("-f", "--size_of_figure", action="store", dest="size_of_figure", type=figsize_from_str,
                    default=(30, 30),
                    help="Size of figure in inches. X and Y values should be separated by comma. Default: 30,30")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=list_from_str,
                    default=["svg", "eps", "pdf", "png", "jpg"],
                    help="Comma-separated list of formats (supported by matplotlib) of output "
                         "figure.Default: svg,eps,pdf,png,jpg")
parser.add_argument("-a", "--annotations", action="store", dest="annotations", required=True,
                    help=".gff file with annotations")
parser.add_argument("-b", "--black_list", action="store", dest="black_list",
                    help="List of annotations to skip")

parser.add_argument("-s", "--synonyms", action="store", dest="synonyms", type=parse_synonyms,
                    help="Set of synonyms. MUST BE i following form: "
                         "<annotation_name_1>:<synonym_to_use>,<annotation_name_2>:<synonym_to_use>")

parser.add_argument("-c", "--combine_mixed_annotations", action="store_true", dest="combine_mixed", default=False,
                    help="Combine mixed annotations")

args = parser.parse_args()

mutations = CollectionVCF(from_file=True, in_file=args.input)
annotations_dict = SeqIO.to_dict(GFF.parse(args.annotations))

mutations.get_location(annotations_dict, use_synonym=True if args.synonyms else False, synonym_dict=args.synonyms)
mutations.location_pie(pie_name="Location of variants", annotation_colors=[],
                       dpi=args.dpi, figsize=args.size_of_figure,
                       explode=True, annotation_black_list=args.black_list,
                       allow_several_counts_of_record=False,
                       pie_prefix=args.output_prefix,
                       full_genome_pie_prefix=args.output_prefix + "_full_genome",
                       counts_filename="location_counts.t",
                       plot_dir="variant_location_pie",
                       counts_dir="location_counts",
                       draw_percents_on_single_pie=False,
                       combine_mixed=args.combine_mixed,
                       extension_list=args.output_formats)

"""
annotation_synonym_dict = {"three_prime_UTR": "3'_UTR",
                           "five_prime_UTR": "5'_UTR",
                           "snoRNA": "ncRNA",
                           "snRNA": "ncRNA"
                           }
annotation_black_list = ["gene", "region", "ARS", "long_terminal_repeat",
                         "noncoding_exon", "intron", "repeat_region", "telomere", "gene_cassette",
                         "five_prime_UTR_intron", "LTR_retrotransposon"]

gene,region,ARS,long_terminal_repeat,noncoding_exon,intron,repeat_region,telomere,gene_cassette,five_prime_UTR_intron,LTR_retrotransposon
three_prime_UTR:3'_UTR,five_prime_UTR:5'_UTR,snoRNA:ncRNA,snRNA:ncRNA
"""






