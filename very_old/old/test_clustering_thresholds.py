#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from MACE.Parsers.VCF import CollectionVCF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with mutations.")
parser.add_argument("-s", "--sample_name", action="store", dest="sample_name", default="unknown_sample",
                    help="Name of sample")
parser.add_argument("-y", "--testing_directory", action="store", dest="test_dir", default="threshold_test",
                    help="Directory where to output results of threshold test")
parser.add_argument("-e", "--extracting_method", action="store", dest="extracting_method", required=True,
                    help="Method used to extract clusters")
parser.add_argument("-p", "--scaffold_prefix", action="store", dest="scaffold_prefix", default="",
                    help="Prefix to write in picture before names of region/scaffold/chromosome. "
                         "Default: no prefix")
parser.add_argument("-d", "--distance_type", action="store", dest="distance_type", default="average",
                    help="Method used to calculate distance between clusters. Default: average")
parser.add_argument("-c", "--count_singletons", action="store_true", dest="count_singletons",
                    help="Draw plot of number of all clusters including singletons. Don't use it for samples "
                         "with low density of mutations")
parser.add_argument("-n", "--min_threshold", action="store", dest="min_threshold", required=True, type=int,
                    help="Minimun threshold for extracting of clusters.")
parser.add_argument("-x", "--max_threshold", action="store", dest="max_threshold", required=True, type=int,
                    help="Maximum threshold for extracting of clusters.")
parser.add_argument("-u", "--number_of_tests", action="store", dest="number_of_tests", required=True, type=int,
                    help="Number of tests")


args = parser.parse_args()

mutations = CollectionVCF(in_file=args.input, from_file=True)

mutations.test_thresholds(extracting_method=args.extracting_method,
                          threshold=(args.min_threshold, args.max_threshold, args.number_of_tests),
                          cluster_distance=args.distance_type,
                          dendrogramm_max_y=2000,
                          sample_name=args.sample_name,
                          save_clustering=False,
                          testing_dir=args.test_dir,
                          count_singletons=args.count_singletons,
                          scaffold_prefix=args.scaffold_prefix)
