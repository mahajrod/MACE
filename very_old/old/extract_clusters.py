#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from MACE.Parsers.VCF import CollectionVCF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with mutations.")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-s", "--sample_name", action="store", dest="sample_name", default="unknown_sample",
                    help="Name of sample")
parser.add_argument("-y", "--clustering_directory", action="store", dest="clust_dir", default="clustering",
                    help="Directory where to output additional data about clustering")
parser.add_argument("-r", "--threshold", action="store", dest="threshold", required=True,
                    help="Threshold for extractig clusters. Depends on extraction method.")
parser.add_argument("-e", "--extracting_method", action="store", dest="extracting_method", required=True,
                    help="Method used to extract clusters")
parser.add_argument("-d", "--distance_type", action="store", dest="distance_type", default="average",
                    help="Method used to calculate distance between clusters. Default: average")

args = parser.parse_args()

mutations = CollectionVCF(in_file=args.input, from_file=True)

clusters = mutations.get_clusters(sample_name=args.sample_name, save_clustering=True,
                                  extracting_method=args.extracting_method,
                                  threshold=args.threshold, cluster_distance=args.distance_type,
                                  dendrogramm_max_y=2000, dendrogramm_color_threshold=1000,
                                  clustering_dir=args.clust_dir)
clusters.write(args.output_prefix + ".ccf")
clusters.write_gff(args.output_prefix + ".gff")
