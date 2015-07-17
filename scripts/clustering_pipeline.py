#!/usr/bin/env python2
import os
import argparse
from collections import OrderedDict

from BCBio import GFF

from MACE.Parsers.VCF import ReferenceGenome, CollectionVCF, ref_alt_variants
from CustomCollections.GeneralCollections import TwoLvlDict


def check_location(feature, pos):
    return True if pos in feature else False


def filter_by_power_05(record):
    return True if record.info_dict['Power'] >= 0.05 else False


def filter_by_power_10(record):
    return True if record.info_dict['Power'] >= 0.10 else False

parser = argparse.ArgumentParser()

parser.add_argument("-r", "--reference_genome", action="store", dest="reference", required=True,
                    help="Fasta file with reference genome")
parser.add_argument("-i", "--reference_genome_index", action="store", dest="reference_index",
                    help="Index of reference genome. Created if absent")
parser.add_argument("-s", "--sample_name", action="store", dest="sample_name", default="unknown_sample",
                    help="Name of sample")
parser.add_argument("-f", "--vcf_file", action="store", dest="vcf_file",
                    help="Vcf file with SNVs")
parser.add_argument("-a", "--annotations", action="store", dest="annotations", required=True,
                    help="Gff file with annotations of reference genome")
parser.add_argument("-m", "--masking", action="store", dest="masking", required=True,
                    help="Gff file with masked regions")
parser.add_argument("-d", "--threshold", action="store", dest="threshold", default=1000, type=int,
                    help="Threshold for extractig clusters. Depends on extraction method.")
parser.add_argument("-y", "--clustering_directory", action="store", dest="clust_dir", default="clustering",
                    help="Directory where to output additional data about clustering")

args = parser.parse_args()


index_file = args.reference_index if args.reference_index else "%s.idx" % (".".join(args.reference.split(".")[:-1]))
reference = ReferenceGenome(args.reference, index_file=index_file)

sample = args.sample_name

clustering_dir = args.clust_dir
distance_threshold = args.threshold
reference.find_gaps()
min_cluster_size = 3

annotations_dict = {}
annotation_synonym_dict = {"three_prime_UTR": "3'_UTR",
                           "five_prime_UTR": "5'_UTR",
                           "snoRNA": "ncRNA",
                           "snRNA": "ncRNA"
                           }
annotation_black_list = ["gene", "region", "ARS", "long_terminal_repeat",
                         "noncoding_exon", "intron", "repeat_region", "telomere", "gene_cassette",
                         "five_prime_UTR_intron"]
with open(args.annotations) as gff_fd:
    for record in GFF.parse(gff_fd):
        annotations_dict[record.id] = record

bad_region_dict = {}
with open(args.masking) as gff_fd:
    for record in GFF.parse(gff_fd):
        bad_region_dict[record.id] = record

statistics_dict = TwoLvlDict(OrderedDict({}))

print("Handling %s" % sample)
statistics_dict[sample] = OrderedDict({})

os.system("mkdir -p %s" % clustering_dir)

mutations = CollectionVCF(in_file=args.vcf_file if args.vcf_file else "%s.vcf" % args.sample_name,
                          from_file=True)

mutations.get_location(annotations_dict, use_synonym=True, synonym_dict=annotation_synonym_dict)
mutations.set_location_flag(bad_region_dict, check_location, "BR")
mutations.check_by_ref_and_alt(ref_alt_variants["deaminases"], "DA", description="Deaminase-like variant")

raw_mutations_counts = len(mutations)
print("Totaly %i mutations" % raw_mutations_counts)
statistics_dict[sample]["raw"] = raw_mutations_counts

sample_adjusted = sample + "_adjusted"
clusters = mutations.get_clusters(sample_name=sample, save_clustering=True,
                                  extracting_method="distance",
                                  threshold=distance_threshold, cluster_distance='average',
                                  dendrogramm_max_y=2000, dendrogramm_color_threshold=1000,
                                  clustering_dir=clustering_dir)
print("Clustering...")
clusters.subclustering()
clusters.statistics(filename="%s/%s_cluster_size_distribution.svg" % (clustering_dir, sample))
clusters.write("%s/%s_raw.ccf" % (clustering_dir, sample))

clusters.adjust(border_limit=None, min_size_to_adjust=1, remove_border_subclusters=True, remove_size_limit=1)
clusters.statistics(filename="%s/%s_cluster_size_distribution.svg" % (clustering_dir, sample_adjusted))

print("Filtering...")
clusters.check_location()
if "HAP" not in sample:
    clusters.check_flags(["DA"], mismatch_list=[1], expression_list=["record.count_samples() <= 1"],
                         remove_mismatch_list=[True])
clusters.get_location(annotations_dict, use_synonym=True, synonym_dict=annotation_synonym_dict)
clusters.write("%s/%s.ccf" % (clustering_dir, sample_adjusted))

filtered_clusters, filtered_out_clusters = clusters.filter_by_flags(black_flag_list=["IP", "BR"])
filtered_clusters.write("%s/%s_not_in_br_no_id.ccf" % (clustering_dir, sample_adjusted))

filtered_out_clusters.write("%s/%s_in_br_id.ccf" % (clustering_dir, sample_adjusted))
filtered_clusters.statistics(filename="%s/%s_not_in_br_no_id_cluster_size_distribution.svg" % (clustering_dir, sample_adjusted))
filtered_out_clusters.statistics(filename="%s/%s_in_br_id_cluster_size_distribution.svg" % (clustering_dir, sample_adjusted))

if "HAP" not in sample:
    filtered_clusters, filtered_out_clusters = filtered_clusters.filter_by_flags(white_flag_list=["DA"])
    filtered_clusters.write("%s/%s_not_in_br_no_id_da.ccf" % (clustering_dir, sample_adjusted))

    filtered_clusters.check_strandness()
    filtered_out_clusters.write("%s/%s_not_in_br_no_id_non_da.ccf" % (clustering_dir, sample_adjusted))


    filtered_clusters.statistics(filename="%s/%s_not_in_br_no_id_da_cluster_size_distribution.svg" % (clustering_dir,sample_adjusted))
    filtered_out_clusters.statistics(filename="%s/%s_not_in_br_no_id_non_da_distribution.svg" % (clustering_dir, sample_adjusted))

    filtered_clusters.heatmap_statistics(filename="%s/%s_not_in_br_no_id_da_heatmap_statistics.svg" % (clustering_dir,sample_adjusted),
                                         additional_data=("Median", "Mean", "Power", "Homogeneity"))
else:
    filtered_clusters.heatmap_statistics(filename="%s/%s_not_in_br_no_id_heatmap_statistics.svg" % (clustering_dir,sample_adjusted),
                                         additional_data=("Median", "Mean", "Power"))


cluster_mutations = filtered_clusters.extract_vcf()
cluster_mutations.write("%s/%s_cluster_mutations.vcf" % (clustering_dir, sample_adjusted))
statistics_dict[sample]["cluster_mutations"] = len(cluster_mutations)
filtered_clusters, filtered_out_clusters = filtered_clusters.filter_by_size(min_size=min_cluster_size)
filtered_clusters.write("%s/%s_size_3+.ccf" % (clustering_dir, sample_adjusted))

filtered_out_clusters.write("%s/%s_size_less_3.ccf" % (clustering_dir, sample_adjusted))
"""
if "HAP" not in sample:
    filtered_clusters.heatmap_statistics(filename="%s/%s_3+_not_in_br_no_id_da_heatmap_statistics.svg" % (clustering_dir,sample_adjusted),
                                         additional_data=("Median", "Mean", "Power", "Homogeneity"))
    #filtered_out_clusters.heatmap_statistics(filename="%s/%s_less_3_not_in_br_no_id_non_da_heatmap_statistics.svg" % (clustering_dir, sample_adjusted),
    #                                         additional_data=("Median", "Mean", "Power", "Homogeneity"))
else:
    filtered_clusters.heatmap_statistics(filename="%s/%s_3+_not_in_br_no_id_heatmap_statistics.svg" % (clustering_dir,sample_adjusted),
                                         additional_data=("Median", "Mean", "Power"))
    #filtered_out_clusters.heatmap_statistics(filename="%s/%s_less_3_not_in_br_no_id_heatmap_statistics.svg" % (clustering_dir, sample_adjusted),
    #                                         additional_data=("Median", "Mean", "Power"))
"""
cluster_mutations = filtered_clusters.extract_vcf()
cluster_mutations.write("%s/%s_3+_cluster_mutations.vcf" % (clustering_dir, sample_adjusted))
statistics_dict[sample]["cluster_3+_mutations"] = len(cluster_mutations)

filtered, filtered_out = filtered_clusters.filter(filter_by_power_05)
filtered.write("%s/%s_adjusted_size_3+_power_0.05+.ccf" % (clustering_dir, sample))
filtered_out.write("%s/%s_adjusted_size_3+_power_less_0.05.ccf" % (clustering_dir, sample))
"""
if "HAP" not in sample:
    filtered.heatmap_statistics(filename="%s/%s_3+_power_0.05+_heatmap_statistics.svg" % (clustering_dir,sample_adjusted),
                                         additional_data=("Median", "Mean", "Power", "Homogeneity"))
    #filtered_out.heatmap_statistics(filename="%s/%s_3+_power_less_0.05_heatmap_statistics.svg" % (clustering_dir, sample_adjusted),
    #                                         additional_data=("Median", "Mean", "Power", "Homogeneity"))
else:
    filtered.heatmap_statistics(filename="%s/%s_3+_power_0.05+_heatmap_statistics.svg" % (clustering_dir,sample_adjusted),
                                         additional_data=("Median", "Mean", "Power"))
    #filtered_out.heatmap_statistics(filename="%s/%s_3+_power_less_0.05_heatmap_statistics.svg" % (clustering_dir, sample_adjusted),
    #                                         additional_data=("Median", "Mean", "Power"))
"""
filtered, filtered_out = filtered.filter(filter_by_power_10)
filtered.write("%s/%s_adjusted_size_3+_power_0.1+.ccf" % (clustering_dir, sample))
filtered_out.write("%s/%s_adjusted_size_3+_power_0.05+_less_0.1.ccf" % (clustering_dir, sample))
"""
if "HAP" not in sample:
    filtered.heatmap_statistics(filename="%s/%s_3+_power_0.10+_heatmap_statistics.svg" % (clustering_dir,sample_adjusted),
                                         additional_data=("Median", "Mean", "Power", "Homogeneity"))
    #filtered_out.heatmap_statistics(filename="%s/%s_3+_power_0.05+_less_0.1_heatmap_statistics.svg" % (clustering_dir, sample_adjusted),
    #                                         additional_data=("Median", "Mean", "Power", "Homogeneity"))
else:
    filtered.heatmap_statistics(filename="%s/%s_3+_power_0.10+_heatmap_statistics.svg" % (clustering_dir,sample_adjusted),
                                         additional_data=("Median", "Mean", "Power"))
    #filtered_out.heatmap_statistics(filename="%s/%s_3+_power_0.05+_less_0.1_heatmap_statistics.svg" % (clustering_dir, sample_adjusted),
    #                                         additional_data=("Median", "Mean", "Power"))
"""
statistics_dict.write(out_filename="%s/%s_mutation_count_statistics.t" % (clustering_dir, sample))