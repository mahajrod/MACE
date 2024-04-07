#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os, sys
import argparse

import numpy as np

from MACE.General.File import split_filename, make_list_of_path_to_files
from MACE.Parsers.VCF import CollectionVCF


def vcf_filter(filename):
    return True if filename[-4:] == ".vcf" else False


def is_homozygous(record):
    return record.is_homozygous()

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output", action="store", dest="output", default="FS_distributions",
                    help="Directory to write output files")

parser.add_argument("-i", "--input_vcf", action="store", dest="input_vcf", type=lambda s: s.split(","),
                    help="Comma-separated list of vcf files or directories containing them",  required=True)

parser.add_argument("-e", "--extension_list", action="store", dest="extension_list", type=lambda s: s.split(","),
                    default=[".png"],
                    help="Comma-separated list of extensions of figures. Default: .png")
args = parser.parse_args()

files_list = sorted(make_list_of_path_to_files(args.input_vcf, vcf_filter))

try:
    os.mkdir(args.output)
except OSError:
    pass

bins = np.arange(0, 66, 5)

for filename in files_list:
    if args.output != "stdout":
        print("Drawing distribution of FS in %s ..." % filename)
    directory, prefix, extension = split_filename(filename)
    variants = CollectionVCF(from_file=True, in_file=filename)

    variants.draw_info_distribution("FS", is_homozygous, outfile_prefix="%s/%s" % (args.output, prefix),
                                    extension_list=args.extension_list, bins=bins)

