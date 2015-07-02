#!/usr/bin/env python2
import os
from collections import OrderedDict

from BCBio import GFF

from MACE.Parsers.VCF import ReferenceGenome, CollectionVCF, ref_alt_variants
from MACE.Parsers.CCF import CollectionCCF
from Parsers.GFF import CollectionGFF
from CustomCollections.GeneralCollections import TwoLvlDict


def check_location(feature, pos):
    return True if pos in feature else False


def filter_by_power_05(record):
    return True if record.info_dict['Power'] >= 0.05 else False


def filter_by_power_10(record):
    return True if record.info_dict['Power'] >= 0.10 else False

if __name__ == "__main__":

    workdir = "analyse/"
    try:
        os.mkdir(workdir)
    except:
        pass
    reference = ReferenceGenome("LAN210_v0.10m.fasta",
                                index_file="LAN210_v0.10m.idx")

    sample_set_names_list = ["PmCDA1_3d"
                             ]

    clustering_dir = "clustering"
    rainfall_dir = "rainfall"
    distance_threshold = 1000
    reference.find_gaps()
    min_cluster_size = 3

    bad_regions_file = "LAN210_v0.10m_masked_all_not_in_good_genes.gff"
    """
    bad_regions = CollectionGFF(input_file=bad_regions_file,
                                from_file=True)
    """
    gff_file = "merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"
    annotations_dict = {}
    annotation_synonym_dict = {"three_prime_UTR": "3'_UTR",
                               "five_prime_UTR": "5'_UTR",
                               "snoRNA": "ncRNA",
                               "snRNA": "ncRNA"
                               }
    annotation_black_list = ["gene", "region", "ARS", "long_terminal_repeat",
                             "noncoding_exon", "intron", "repeat_region", "telomere", "gene_cassette",
                             "five_prime_UTR_intron"]
    with open(gff_file) as gff_fd:
        for record in GFF.parse(gff_fd):
            annotations_dict[record.id] = record

    bad_region_dict = {}
    with open(bad_regions_file) as gff_fd:
        for record in GFF.parse(gff_fd):
            bad_region_dict[record.id] = record

    statistics_dict = TwoLvlDict(OrderedDict({}))
    for sample_set_name in sample_set_names_list:
        print("Handling %s" % sample_set_name)
        statistics_dict[sample_set_name] = OrderedDict({})
        os.chdir(workdir)
        os.system("mkdir -p %s" % sample_set_name)
        os.chdir(sample_set_name)
        os.system("mkdir -p %s" % clustering_dir)

        mutations = CollectionVCF(in_file="../../%s_SNP.vcf" % sample_set_name,
                                  from_file=True)

        mutations.get_location(annotations_dict, use_synonym=True, synonym_dict=annotation_synonym_dict)
        mutations.set_location_flag(bad_region_dict, check_location, "BR")
        mutations.check_by_ref_and_alt(ref_alt_variants["deaminases"], "DA", description="Deaminase-like variant")

        raw_mutations_counts = len(mutations)
        print("Totaly %i mutations" % raw_mutations_counts)
        statistics_dict[sample_set_name]["raw"] = raw_mutations_counts

        sample_set_name_adjusted = sample_set_name + "_adjusted"
        clusters = mutations.get_clusters(sample_name=sample_set_name, save_clustering=True,
                                          extracting_method="distance",
                                          threshold=distance_threshold, cluster_distance='average',
                                          dendrogramm_max_y=2000, dendrogramm_color_threshold=1000,
                                          clustering_dir=clustering_dir)
        clusters.subclustering()
        clusters.statistics(filename="%s/%s_cluster_size_distribution.svg" % (clustering_dir, sample_set_name))
        clusters.write("%s/%s_raw.ccf" % (clustering_dir, sample_set_name))

        clusters.adjust(border_limit=None, min_size_to_adjust=1, remove_border_subclusters=True, remove_size_limit=1)
        clusters.statistics(filename="%s/%s_cluster_size_distribution.svg" % (clustering_dir, sample_set_name_adjusted))

        clusters.check_location()
        if "HAP" not in sample_set_name:
            clusters.check_flags(["DA"], mismatch_list=[1], expression_list=["record.count_samples() <= 1"],
                                 remove_mismatch_list=[True])
        clusters.get_location(annotations_dict, use_synonym=True, synonym_dict=annotation_synonym_dict)
        clusters.write("%s/%s.ccf" % (clustering_dir, sample_set_name_adjusted))

        filtered_clusters, filtered_out_clusters = clusters.filter_by_flags(black_flag_list=["IP", "BR"])
        filtered_clusters.write("%s/%s_not_in_br_no_id.ccf" % (clustering_dir, sample_set_name_adjusted))

        filtered_out_clusters.write("%s/%s_in_br_id.ccf" % (clustering_dir, sample_set_name_adjusted))
        filtered_clusters.statistics(filename="%s/%s_not_in_br_no_id_cluster_size_distribution.svg" % (clustering_dir, sample_set_name_adjusted))
        filtered_out_clusters.statistics(filename="%s/%s_in_br_id_cluster_size_distribution.svg" % (clustering_dir, sample_set_name_adjusted))

        if "HAP" not in sample_set_name:
            filtered_clusters, filtered_out_clusters = filtered_clusters.filter_by_flags(white_flag_list=["DA"])
            filtered_clusters.write("%s/%s_not_in_br_no_id_da.ccf" % (clustering_dir, sample_set_name_adjusted))

            filtered_clusters.check_strandness()
            filtered_out_clusters.write("%s/%s_not_in_br_no_id_non_da.ccf" % (clustering_dir, sample_set_name_adjusted))


            filtered_clusters.statistics(filename="%s/%s_not_in_br_no_id_da_cluster_size_distribution.svg" % (clustering_dir,sample_set_name_adjusted))
            filtered_out_clusters.statistics(filename="%s/%s_not_in_br_no_id_non_da_distribution.svg" % (clustering_dir, sample_set_name_adjusted))

            filtered_clusters.heatmap_statistics(filename="%s/%s_not_in_br_no_id_da_heatmap_statistics.svg" % (clustering_dir,sample_set_name_adjusted),
                                                 additional_data=("Median", "Mean", "Power", "Homogeneity"))
        else:
            filtered_clusters.heatmap_statistics(filename="%s/%s_not_in_br_no_id_heatmap_statistics.svg" % (clustering_dir,sample_set_name_adjusted),
                                                 additional_data=("Median", "Mean", "Power"))


        cluster_mutations = filtered_clusters.extract_vcf()
        cluster_mutations.write("%s/%s_cluster_mutations.vcf" % (clustering_dir, sample_set_name_adjusted))
        statistics_dict[sample_set_name]["cluster_mutations"] = len(cluster_mutations)
        filtered_clusters, filtered_out_clusters = filtered_clusters.filter_by_size(min_size=min_cluster_size)
        filtered_clusters.write("%s/%s_size_3+.ccf" % (clustering_dir, sample_set_name_adjusted))

        filtered_out_clusters.write("%s/%s_size_less_3.ccf" % (clustering_dir, sample_set_name_adjusted))

        if "HAP" not in sample_set_name:
            filtered_clusters.heatmap_statistics(filename="%s/%s_3+_not_in_br_no_id_da_heatmap_statistics.svg" % (clustering_dir,sample_set_name_adjusted),
                                                 additional_data=("Median", "Mean", "Power", "Homogeneity"))
            filtered_out_clusters.heatmap_statistics(filename="%s/%s_less_3_not_in_br_no_id_non_da_heatmap_statistics.svg" % (clustering_dir, sample_set_name_adjusted),
                                                     additional_data=("Median", "Mean", "Power", "Homogeneity"))
        else:
            filtered_clusters.heatmap_statistics(filename="%s/%s_3+_not_in_br_no_id_heatmap_statistics.svg" % (clustering_dir,sample_set_name_adjusted),
                                                 additional_data=("Median", "Mean", "Power"))
            filtered_out_clusters.heatmap_statistics(filename="%s/%s_less_3_not_in_br_no_id_heatmap_statistics.svg" % (clustering_dir, sample_set_name_adjusted),
                                                     additional_data=("Median", "Mean", "Power"))

        cluster_mutations = filtered_clusters.extract_vcf()
        cluster_mutations.write("%s/%s_3+_cluster_mutations.vcf" % (clustering_dir, sample_set_name_adjusted))
        statistics_dict[sample_set_name]["cluster_3+_mutations"] = len(cluster_mutations)

        filtered, filtered_out = filtered_clusters.filter(filter_by_power_05)
        filtered.write("%s/%s_adjusted_size_3+_power_0.05+.ccf" % (clustering_dir, sample_set_name))
        filtered_out.write("%s/%s_adjusted_size_3+_power_less_0.05.ccf" % (clustering_dir, sample_set_name))

        if "HAP" not in sample_set_name:
            filtered.heatmap_statistics(filename="%s/%s_3+_power_0.05+_heatmap_statistics.svg" % (clustering_dir,sample_set_name_adjusted),
                                                 additional_data=("Median", "Mean", "Power", "Homogeneity"))
            #filtered_out.heatmap_statistics(filename="%s/%s_3+_power_less_0.05_heatmap_statistics.svg" % (clustering_dir, sample_set_name_adjusted),
            #                                         additional_data=("Median", "Mean", "Power", "Homogeneity"))
        else:
            filtered.heatmap_statistics(filename="%s/%s_3+_power_0.05+_heatmap_statistics.svg" % (clustering_dir,sample_set_name_adjusted),
                                                 additional_data=("Median", "Mean", "Power"))
            #filtered_out.heatmap_statistics(filename="%s/%s_3+_power_less_0.05_heatmap_statistics.svg" % (clustering_dir, sample_set_name_adjusted),
            #                                         additional_data=("Median", "Mean", "Power"))

        filtered, filtered_out = filtered.filter(filter_by_power_10)
        filtered.write("%s/%s_adjusted_size_3+_power_0.1+.ccf" % (clustering_dir, sample_set_name))
        filtered_out.write("%s/%s_adjusted_size_3+_power_0.05+_less_0.1.ccf" % (clustering_dir, sample_set_name))

        if "HAP" not in sample_set_name:
            filtered.heatmap_statistics(filename="%s/%s_3+_power_0.10+_heatmap_statistics.svg" % (clustering_dir,sample_set_name_adjusted),
                                                 additional_data=("Median", "Mean", "Power", "Homogeneity"))
            #filtered_out.heatmap_statistics(filename="%s/%s_3+_power_0.05+_less_0.1_heatmap_statistics.svg" % (clustering_dir, sample_set_name_adjusted),
            #                                         additional_data=("Median", "Mean", "Power", "Homogeneity"))
        else:
            filtered.heatmap_statistics(filename="%s/%s_3+_power_0.10+_heatmap_statistics.svg" % (clustering_dir,sample_set_name_adjusted),
                                                 additional_data=("Median", "Mean", "Power"))
            #filtered_out.heatmap_statistics(filename="%s/%s_3+_power_0.05+_less_0.1_heatmap_statistics.svg" % (clustering_dir, sample_set_name_adjusted),
            #                                         additional_data=("Median", "Mean", "Power"))

    #statistics_dict.write(out_filename=workdir + "mutation_count_statistics.t")