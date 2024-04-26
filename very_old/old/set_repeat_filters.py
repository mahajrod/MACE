#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Bio import SeqIO
from BCBio import GFF

from MACE.Parsers.VCF import CollectionVCF


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_vcf", action="store", dest="input_vcf", required=True,
                    help="Input vcf file")
parser.add_argument("-o", "--output_vcf", action="store", dest="output_vcf", required=True,
                    help="Output vcf file")
parser.add_argument("-r", "--repeat_gff", action="store", dest="repeat_gff", required=True,
                    help="Gff with repeats")
parser.add_argument("-f", "--filter_name", action="store", dest="filter_name", default="in_repeat",
                    help="Name of filter to be set. Default - 'in_repeat'")
args = parser.parse_args()

annotation_dict = SeqIO.to_dict(GFF.parse(args.repeat_gff))

vcf_collection = CollectionVCF(from_file=True, in_file=args.input_vcf)
vcf_collection.set_filter_by_intersection_with_feature(annotation_dict, args.filter_name, mode="cross",
                                                       feature_type_black_list=[])
vcf_collection.write(args.output_vcf)


