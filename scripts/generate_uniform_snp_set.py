#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Bio import SeqIO
from BCBio import GFF

from MACE.Parsers.VCF import ReferenceGenome


def parse_substititions(s):
    if not s:
        return None
    temp = s.split("-")
    substitution_dict = {}
    for entry in temp:
        t = entry.split(":")
        if len(t) > 1:
            substitution_dict[t[0]] = t[1].split(",")
        else:
            substitution_dict[t[0]] = None
    return substitution_dict



parser = argparse.ArgumentParser()

parser.add_argument("-r", "--reference_genome", action="store", dest="reference_genome", required=True,
                    help="File with reference genome")
parser.add_argument("-g", "--masking_gff", action="store", dest="masking_gff",
                    help=".gff file with masked regions")
parser.add_argument("-i", "--reference_genome_index", action="store", dest="ref_gen_idx", default="refgen.idx",
                    help="Index of reference genome")
parser.add_argument("-m", "--parsing_mode", action="store", dest="parsing_mode", default="index_db",
                    help="Parsing mode of reference genome. Allowed: to_dict, index, index_db. The fastest is to_dict, "
                         "the lowest memory consuming - index_db. Default: index_db")
parser.add_argument("-b", "--region_black_list", action="store", dest="black_list", type=lambda s: s.split(","),
                    default=[],
                    help="Comma-separated ist of region names in genome to be not mutated")
parser.add_argument("-v", "--out_vcf", action="store", dest="out_vcf", required=True,
                    help=".vcf with snp set")
parser.add_argument("-n", "--number_of_mutations", action="store", dest="mut_number", type=int, default=10000,
                    help="Number of mutations in set")
parser.add_argument("-z", "--zygoty", action="store", dest="zygoty", default="homo",
                    help="Zygoty of mutations in set. At moment only 100%% heterozygous or 100%% "
                         "homozygous sets can be generated. "
                         "Allowed values: homo, hetero. Default: homo")
parser.add_argument("-s", "--substitutions", action="store", dest="substitutions", type=parse_substititions,
                    help="Set of substitution. MUST BE i following form: "
                         "<ref_base_1>:<comma-separetad_alternatives>-<ref_base_2>:<comma-separetad_alternatives> "
                         "Alternatives can be not set. If so all possible variants will be choosen as "
                         "a set of alternatives."
                         "If no reference bases was set - sites with all four bases will be considered as mutation sites. "
                         "Example:  G:T,C,A-T-A:G ")

args = parser.parse_args()

masked_regions_dict = SeqIO.to_dict(GFF.parse(args.masking_gff))

reference_genome = ReferenceGenome(args.reference_genome, masked_regions=masked_regions_dict,
                                   index_file=args.ref_gen_idx, filetype="fasta", mode=args.parsing_mode,
                                   black_list=args.black_list)
reference_genome.generate_snp_set(args.mut_number, substitution_dict=args.substitutions,
                                  zygoty=args.zygoty, out_vcf=args.out_vcf)

