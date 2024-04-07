#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from Bio import SeqIO

from MACE.Parsers.VCF import CollectionVCF


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_vcf", action="store", dest="input_vcf", required=True,
                    help="Input vcf file")
parser.add_argument("-o", "--output_vcf", action="store", dest="output_vcf", required=True,
                    help="Output vcf file")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Reference file")
parser.add_argument("-f", "--filter_name", action="store", dest="filter_name", default="in_repeat",
                    help="Name of filter to be set. Default - 'indel_in_homopolymer'")
parser.add_argument("-m", "--min_homopolymer_len", action="store", dest="min_homopolymer_len",
                    default=4, type=int,
                    help="Minimum length of homopolymer to set a filter. Default - 4")
args = parser.parse_args()

reference_dict = SeqIO.index_db("tmp.idx", args.reference, format="fasta")
vcf_collection = CollectionVCF(from_file=True, in_file=args.input_vcf)
vcf_collection.set_filter_for_indels_in_homopolymers(reference_dict, min_homopolymer_len=args.min_homopolymer_len,
                                                     filter_name=args.filter_name)
vcf_collection.write(args.output_vcf)
os.remove("tmp.idx")


