#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from collections import OrderedDict
from BCBio import GFF
from MACE.Parsers.VCF import CollectionVCF, ReferenceGenome
from MACE.Routines import DrawingRoutines
from MACE.General.GeneralCollections import IdList

parser = argparse.ArgumentParser()

parser.add_argument("-a", "--vcf_a", action="store", dest="vcf_a", required=True,
                    help="First vcf file")
parser.add_argument("-b", "--vcf_b", action="store", dest="vcf_b", required=True,
                    help="Second vcf file")
parser.add_argument("--name_a", action="store", dest="name_a", default="vcf_a",
                    help="Sample name for first vcf file. Default: vcf_a")
parser.add_argument("--name_b", action="store", dest="name_b", default="vcf_b",
                    help="Sample name for second vcf file. Default: vcf_b")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-r", "--reference_genome", action="store", dest="reference",
                    help="Fasta file with reference genome, required to draw gaps and chromosomes")

parser.add_argument("-w", "--window_size", action="store", dest="window_size", default=100000, type=int,
                    help="Size of the windows Default: 100000")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", default=None, type=int,
                    help="Step of the sliding windows. Default: window size, i.e windows are staking")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="index_db",
                    help="Parsing mode for input sequence file. "
                         "Possible variants: 'index_db'(default), 'index', 'parse'")
parser.add_argument("-x", "--scaffold_white_list", action="store", dest="scaffold_white_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of the only scaffolds to draw. Default: all")
parser.add_argument("-u", "--scaffold_black_list", action="store", dest="scaffold_black_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of scaffolds to skip at drawing. Default: not set")
parser.add_argument("-y", "--sort_scaffolds", action="store_true", dest="sort_scaffolds", default=False,
                    help="Order  scaffolds according to their names. Default: False")
parser.add_argument("-z", "--scaffold_ordered_list", action="store", dest="scaffold_ordered_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of scaffolds to draw first and exactly in same order. "
                         "Scaffolds absent in this list are drawn last and in order according to vcf file . "
                         "Default: not set")
parser.add_argument("-m", "--minimal_ratio", action="store", dest="minimal_ratio", default=2, type=float,
                    help="Minimal ratio between densities in regions. Default: 2")



args = parser.parse_args()

count_dict = OrderedDict()

reference = ReferenceGenome(args.reference,
                            masked_regions=None,
                            index_file="refgen.idx",
                            filetype="fasta",
                            mode=args.parsing_mode,
                            black_list=[])

reference.find_gaps(min_gap_length=10)

for sample_name, vcf_file in ((args.name_a, args.vcf_a), (args.name_b, args.vcf_b)):
    count_dict[sample_name] = CollectionVCF(from_file=True, in_file=vcf_file, parse_only_coordinates=True).count_variants_in_windows(args.window_size,
                                                                                                                                     args.window_size if args.window_step is None else args.window_step,
                                                                                                                                     reference.region_length,
                                                                                                                                     ignore_scaffolds_shorter_than_window=True,
                                                                                                                                     output_prefix="%s.%s" % (args.output_prefix, sample_name),
                                                                                                                                     skip_empty_windows=False)

white_set = set(args.scaffold_white_list)
black_set = set(args.scaffold_black_list)

scaffold_set = set()
for sample in count_dict:
    scaffold_set |= set(count_dict[sample])

if white_set:
    scaffold_set = scaffold_set & white_set

if black_set:
    scaffold_set = scaffold_set - black_set

scaffold_list = list(scaffold_set)

if args.sort_scaffolds:
    scaffold_list.sort()

final_scaffold_list = []

if args.scaffold_ordered_list:
    for entry in scaffold_ordered_list:
        final_scaffold_list.append(entry)
        scaffold_list.remove(entry)
    final_scaffold_list = final_scaffold_list + scaffold_list
else:
    final_scaffold_list = scaffold_list

vcf_a_more_variants_file = "%s.higher_a.tab" % args.output_prefix
vcf_b_more_variants_file = "%s.higher_b.tab" % args.output_prefix

vcf_a_no_variants_file = "%s.a_zero.tab" % args.output_prefix
vcf_b_no_variants_file = "%s.b_zero.tab" % args.output_prefix

vcf_a_absent_scaffolds_id_file = "%s.a_absent_scaffolds.ids" % args.output_prefix
vcf_b_absent_scaffolds_id_file = "%s.b_absent_scaffolds.ids" % args.output_prefix

vcf_density_ratio_file = "%s.a_to_b_density_ratio.tab" % args.output_prefix

vcf_a_more_variants_file_fd = open(vcf_a_more_variants_file, "w")
vcf_b_more_variants_file_fd = open(vcf_a_more_variants_file, "w")
vcf_a_no_variants_file_fd = open(vcf_a_no_variants_file, "w")
vcf_b_no_variants_file_fd = open(vcf_b_no_variants_file, "w")

vcf_density_ratio_fd = open(vcf_density_ratio_file, "w")

vcf_a_absent_scaffolds_id_list = IdList()
vcf_b_absent_scaffolds_id_list = IdList()

vcf_a_more_variants_file_fd.write("#scaffold\tstart\tend\twindow_id(0-based)\t%s(N_of_variants)\t%s(N_of_variants)\tratio\n" % (args.name_a, args.name_b))
vcf_b_more_variants_file_fd.write("#scaffold\tstart\tend\twindow_id(0-based)\t%s(N_of_variants)\t%s(N_of_variants)\tratio\n" % (args.name_a, args.name_b))

vcf_a_no_variants_file_fd.write("#scaffold\tstart\tend\twindow_id(0-based)\n")
vcf_b_no_variants_file_fd.write("#scaffold\tstart\tend\twindow_id(0-based)\n")

vcf_density_ratio_fd.write("#scaffold\t%s\t%s\tratio,a/b\n" % (args.name_a, args.name_b))

for scaffold in final_scaffold_list:
    if (scaffold in count_dict[args.name_a]) and (scaffold in count_dict[args.name_b]):
        for i in range(0, len(count_dict[args.name_a][scaffold])):
            start = i*args.window_step + 1
            stop = start + args.window_size - 1

            if (count_dict[args.name_a][scaffold] > 0) and (count_dict[args.name_b][scaffold] > 0):
                ratio = float(count_dict[args.name_a][scaffold])/float(count_dict[args.name_b][scaffold])
                vcf_density_ratio_fd.write("%s\t%i\t%i\t%.3f\n" % (scaffold,
                                                                   count_dict[args.name_a][scaffold],
                                                                   count_dict[args.name_b][scaffold],
                                                                   ratio))
                if ratio > args.minimal_ratio:
                    vcf_a_more_variants_file.write("%s\t%i\t%i\t%i\t%i\t%i\t%.3f\n" % (scaffold, start, stop, i,
                                                                                       count_dict[args.name_a][scaffold],
                                                                                       count_dict[args.name_b][scaffold],
                                                                                       ratio))
                elif ratio < (1.0/float(args.minimal_ratio)):
                    vcf_b_more_variants_file.write("%s\t%i\t%i\t%i\t%i\t%i\t%.3f\n" % (scaffold, start, stop, i,
                                                                                       count_dict[args.name_a][scaffold],
                                                                                       count_dict[args.name_b][scaffold],
                                                                                       ratio))
                    
            elif count_dict[args.name_a][scaffold] == 0:
                vcf_a_no_variants_file_fd.write("%s\t%i\t%i\t%i" % (scaffold, start, stop, i))
            elif count_dict[args.name_b][scaffold] == 0:
                vcf_b_no_variants_file_fd.write("%s\t%i\t%i\t%i" % (scaffold, start, stop, i))

    else:
        if scaffold not in count_dict[args.name_a]:
            vcf_a_absent_scaffolds_id_list.append(scaffold)
        if scaffold not in count_dict[args.name_b]:
            vcf_b_absent_scaffolds_id_list.append(scaffold)




vcf_a_more_variants_file_fd.close()
vcf_b_more_variants_file_fd.close()
vcf_a_no_variants_file_fd.close()
vcf_b_no_variants_file_fd.close()

vcf_density_ratio_fd.close()

vcf_a_absent_scaffolds_id_list.write(vcf_a_absent_scaffolds_id_file)
vcf_b_absent_scaffolds_id_list.write(vcf_b_absent_scaffolds_id_file)