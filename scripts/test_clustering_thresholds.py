#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from RouToolPa.Parsers.VCF import CollectionVCF
from MACE.Routines import StatsVCF, Visualization
from RouToolPa.Collections.General import IdList, SynDict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with variants")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-d", "--distance", action="store", dest="distance", default="average",
                    help="Method to use for calculating distance between clusters. "
                         "Allowed: average(default), single, complete, centroid, weighted, median, ward")
parser.add_argument("-e", "--method", action="store", dest="method", default="distance",
                    help="Method to use for extraction of clusters. "
                         "Allowed: distance(default), inconsistent, maxclust, monocrit, maxclust_monocrit ")

parser.add_argument("-m", "--min_threshold", action="store", dest="min_threshold", type=float,
                    help="Minimum threshold for extraction of clusters")
parser.add_argument("-x", "--max_threshold", action="store", dest="max_threshold", type=float,
                    help="Maximum threshold for extraction of clusters")
parser.add_argument("-s", "--threshold_step", action="store", dest="threshold_step", type=float,
                    help="Threshold step for extraction of clusters")

parser.add_argument("-a", "--scaffold_white_list", action="store", dest="scaffold_white_list", default=[],
                    type=lambda s: IdList(filename=s) if os.path.exists(s) else s.split(","),
                    help="Comma-separated list of the only scaffolds to draw. Default: all")

parser.add_argument("-b", "--scaffold_black_list", action="store", dest="scaffold_black_list", default=[],
                    type=lambda s: IdList(filename=s) if os.path.exists(s) else s.split(","),
                    help="Comma-separated list of scaffolds to skip at drawing. Default: not set")

parser.add_argument("-z", "--scaffold_ordered_list", action="store", dest="scaffold_ordered_list", default=[],
                    type=lambda s: IdList(filename=s) if os.path.exists(s) else s.split(","),
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

args = parser.parse_args()

syn_dict = SynDict(filename=args.scaffold_syn_file,
                   key_index=args.syn_file_key_column,
                   value_index=args.syn_file_value_column)

variants = CollectionVCF(in_file=args.input, parsing_mode="only_coordinates",
                         scaffold_white_list=args.scaffold_white_list,
                         scaffold_black_list=args.scaffold_black_list,
                         scaffold_syn_dict=syn_dict)

linkage_df = StatsVCF.get_linkage_for_hierarchical_clustering(variants.records, method=args.distance, output=None)

cluster_df = StatsVCF.test_clustering_thresholds(linkage_df,
                                                 extracting_method=args.method,
                                                 threshold_tuple=None,
                                                 min_threshold=args.min_threshold,
                                                 max_threshold=args.max_threshold,
                                                 threshold_number=None,
                                                 threshold_step=args.threshold_step,
                                                 output_prefix=None)

Visualization.plot_clustering_threshold_tests(cluster_df, args.output_prefix,
                                              scaffold_order_list=args.scaffold_ordered_list,
                                              extensions=("png", ), suptitle="Test of clustering thresholds")
