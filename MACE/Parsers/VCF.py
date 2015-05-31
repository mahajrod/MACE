#!/usr/bin/env python
"""
VCF Parser Module
"""
__author__ = 'Sergei F. Kliver'

import os
import re
from math import sqrt
from copy import deepcopy
from collections import OrderedDict

import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, inconsistent, cophenet, fcluster

import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from MACE.Parsers.Abstract import Record, Collection, Metadata, Header


ref_alt_variants = {"desaminases": [("C", ["T"]), ("G", ["A"])]
                    }


class RecordVCF(Record):
    """
    RecordVCF class
    """
    def __init__(self, pos, id, ref, alt_list, qual, filter_list, info_dict, samples_list,
                 description=None, flags=None):
        """

        :param pos:
        :param id:
        :param ref:
        :param alt_list:
        :param qual:
        :param filter_list:
        :param info_dict:
        :param samples_list:
        :param description:
        :param flags:
        :return:
        """

        self.pos = pos                                  #int
        self.id = id                                    #str
        self.ref = ref                                  #str
        self.alt_list = alt_list                        #list, entries are strings
        self.qual = qual                                #real or "."
        self.filter_list = sorted(filter_list)          #list, entries are strings
        self.info_dict = info_dict                      #dict
        self.samples_list = samples_list                #list entries are dicts with keys from format_list and
                                                        #values are lists
        #self.description = description if description else {}
        self.flags = set(flags) if flags is not None else set([])

    def __str__(self):
        """

        :return:
        """
        alt_string = ",".join(self.alt_list)
        filter_string = ";".join(self.filter_list)
        key_string_list = []
        for key in sorted(list(self.info_dict.keys())):
            if self.info_dict[key]:
                key_string_list.append(key + "=" + ",". join(map(lambda x: str(x), self.info_dict[key])))
            else:
                key_string_list.append(key)

        info_string = ";".join(key_string_list)
        info_string = "%s;%s" % (info_string, ";".join(self.flags)) if self.flags else info_string
        for sample in self.samples_list:
            if len(sample.keys()) > 1:
                format_string = ":".join(sample.keys())
                break

        samples_string = "\t".join([":".join([",".join(map(lambda x: str(x), sample[key])) for key in sample.keys()]) for sample in self.samples_list])
        return '\t'.join(map(lambda x: str(x), [self.pos, self.id, self.ref, alt_string,
                                                self.qual, filter_string, info_string, format_string, samples_string]))

    def check_ref_alt_list(self, ref_alt_list, flag):
        """

        :param ref_alt_list:
        :param flag:
        :return:
        """
        # structure of ref_alt_list:  [[ref1,[alt1.1, alt1.M1]], ..., [refN,[altN.1, ..., altN.MN]]]
        self.set_flag(lambda record: (record.ref, record.alt_list) in ref_alt_list, flag)

    def count_samples(self):
        """

        :return:
        """
        #counts samples with variants
        number = 0
        for sample in self.samples_list:
            if sample["GT"][0] != "./." and sample["GT"][0] != "0/0":
                number += 1
        return number

    def set_filter(self, expression, filter_name):
        """
        :param expression:
        :param filter_name:
        :return:
        """
        if expression(self):
            self.filter_list.append(filter_name)
            self.filter_list.sort()

    def add_info(self, info_name, info_value=None):
        """

        :param info_name:
        :param info_value:
        :return:
        """
        value = info_value if isinstance(info_value, list) else [] if info_value is None else [info_value]
        if info_name in self.info_dict:
            self.info_dict[info_name] += value
        else:
            self.info_dict[info_name] = value

    def check_indel(self):
        """
        Checks if record is indel
        """
        if len(self.ref) > 1 or len("".join(self.alt_list)) > len(self.alt_list):
            return True
        return False

    def find_location(self, scaffold, annotation_dict, key="Ftype", strand_key="Fstrand", genes_key="Genes",
                      genes_strand_key="Gstrand", feature_type_black_list=[],
                      use_synonym=False, synonym_dict=None, add_intergenic_label=True):
        """

        :param scaffold:
        :param annotation_dict:
        :param key:
        :param strand_key:
        :param genes_key:
        :param genes_strand_key:
        :param feature_type_black_list:
        :param use_synonym:
        :param synonym_dict:
        :param add_intergenic_label:
        :return:
        """
        """
        This method is written for old variant (with sub_feature)s rather then new (with CompoundLocation)
        id of one SeqRecord in record_dict must be equal to record.pos
        locations will be written to description dictionary of record using "key" as key
        """
        if key not in self.info_dict:
            self.info_dict[key] = set([])

        # strandness values:
        # N - not defined
        # B - both
        # P - plus
        # M - minus
        strands = ["B", "P", "M"]
        if strand_key not in self.info_dict:
            self.info_dict[strand_key] = ["N"]
        for flag_key in (genes_key, genes_strand_key):
            if flag_key not in self.info_dict:
                self.info_dict[flag_key] = []

        for feature in annotation_dict[scaffold].features:
            if feature.type in feature_type_black_list:
                continue

            if (self.pos - 1) in feature:
                if feature.type == "gene" or feature.type == "ncRNA":
                    self.info_dict[genes_key].append(feature.qualifiers["Name"][0])
                    self.info_dict[genes_strand_key].append(strands[feature.strand])

                self.info_dict[key].add(self.get_synonym(feature.type, use_synonym=use_synonym,
                                                         synonym_dict=synonym_dict))
                if self.info_dict[strand_key][0] == "N":
                    self.info_dict[strand_key][0] = strands[feature.strand]
                elif strands[feature.strand] != self.info_dict[strand_key][0]:
                    self.info_dict[strand_key][0] = "B"
            else:
                continue

            for sub_feature in feature.sub_features:
                if sub_feature.type in feature_type_black_list:
                    continue
                if (self.pos - 1) in sub_feature:
                    self.info_dict[key].add(self.get_synonym(sub_feature.type, use_synonym=use_synonym,
                                                             synonym_dict=synonym_dict))
                    if self.info_dict[strand_key][0] == "N":
                        self.info_dict[strand_key][0] = strands[sub_feature.strand]
                    elif strands[sub_feature.strand] != self.info_dict[strand_key][0]:
                        self.info_dict[strand_key][0] = "B"

        if not self.info_dict[genes_key]:
            self.info_dict.pop(genes_key)
            self.info_dict.pop(genes_strand_key)

        if add_intergenic_label and (not self.info_dict[key]): # or ("gene" not in self.info_dict[key])):
            # igc == intergenic
            self.info_dict[key].add("igc")


class MetadataVCF(OrderedDict, Metadata):
    """
    MetadataVCF class
    """
    @staticmethod
    def _split_by_equal_sign(string):
        try:
            index = string.index("=")
        except ValueError:
            #if "=" is not present in string (case of flag type in INFO field)
            return string, None
        return string[:index], string[index+1:]

    @staticmethod
    def _split_by_comma_sign(string):
        index_list = [-1]
        i = 1
        while (i < len(string)):
            if string[i] == "\"":
                i += 1
                while string[i] != "\"":
                    i += 1
            if string[i] == ",":
                index_list.append(i)
            i += 1
        index_list.append(len(string))
        return [string[index_list[j] + 1: index_list[j + 1]] for j in range(0, len(index_list) - 1)]

    def add_metadata(self, line):
        key, value = self._split_by_equal_sign(line[2:].strip())
        if value[0] == "<" and value[-1] == ">":
            value = self._split_by_comma_sign(value[1:-1])
            value_id = self._split_by_equal_sign(value[0])[1]
            value = OrderedDict(self._split_by_equal_sign(entry) for entry in value[1:])
            if key not in self:
                self[key] = OrderedDict({})
            self[key][value_id] = value
        else:
            self[key] = value

    def add_metadata_from_values(self, name, number, ntype, description):
        self[name] = OrderedDict({})
        self[name]["Number"] = number
        self[name]["Type"] = ntype
        self[name]["Description"] = description

    def __str__(self):
        metadata_string = ""
        for key in self:
            if not isinstance(self[key], dict):
                metadata_string += "##%s=%s\n" % (key, self[key])
            else:
                prefix = "##%s=<" % key
                suffix = ">\n"
                for att_id in self[key]:
                    middle = "ID=%s," % att_id + ",".join(["%s=%s" % (param, self[key][att_id][param])
                                                           for param in self[key][att_id]])
                    metadata_string += prefix + middle + suffix
        return metadata_string[:-1]


class HeaderVCF(list, Header):
    """
    HeaderVCF class
    """
    def __str__(self):
        return "#" + "\t".join(self)


class CollectionVCF(Collection):
    """
    CollectionVCF class

    """

    def __init__(self, metadata=None, records_dict=None, header=None, in_file=None, samples=None,
                 from_file=True, external_metadata=None, threads=1):
        """

        :param metadata:
        :param records_dict:
        :param header:
        :param in_file:
        :param samples:
        :param from_file:
        :param external_metadata:
        :param threads:
        :return:
        """
        self.linkage_dict = None
        if from_file:
            self.read(in_file, external_metadata=external_metadata)
        else:
            self.metadata = metadata
            self.records = {} if records_dict is None else records_dict
            self.header = header
            self.samples = samples
        self.scaffold_list = self.scaffolds()
        self.scaffold_length = self.scaffold_len()
        self.number_of_scaffolds = len(self.scaffold_list)
        self.record_index = self.rec_index()
        self.threads = threads

    def read(self, in_file, external_metadata=None):
        """

        :param in_file:
        :param external_metadata:
        :return:
        """
        self.metadata = MetadataVCF()
        self.records = OrderedDict({})
        with open(in_file, "r") as fd:
            for line in fd:
                if line[:2] != "##":
                    self.header = HeaderVCF(line[1:].strip().split("\t"))   # line[1:].strip().split("\t")
                    self.samples = self.header[9:]
                    break
                self.metadata.add_metadata(line)
            for line in fd:
                scaffold, record = self.add_record(line, external_metadata=external_metadata)
                if scaffold not in self.records:
                    self.records[scaffold] = []
                self.records[scaffold].append(record)

    def add_record(self, line, external_metadata=None):
        """

        :param line:
        :param external_metadata:
        :return:
        """
        line_list = line.strip().split("\t")
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_1
        position = int(line_list[1])
        quality = "."
        if quality != line_list[5]:
            quality = float(line_list[5])
        alt_list = line_list[4].split(",")
        filter_list = line_list[6].split(",")          # list, entries are strings

        info_dict = OrderedDict()
        if line_list[7] != ".":
            info_tuple_list = [self._split_by_equal_sign(entry) for entry in line_list[7].split(";")]

            flag_set = set([])
            metadata = self.metadata if self.metadata else external_metadata
            for entry in info_tuple_list:
                # if self.metadata:
                if metadata["INFO"][entry[0]]["Type"] == "Flag":
                    flag_set.add(entry[0]) #info_dict[entry[0]] = []
                elif metadata["INFO"][entry[0]]["Type"] == "Integer":
                    info_dict[entry[0]] = list(map(lambda x: int(x), entry[1].split(",")))
                elif metadata["INFO"][entry[0]]["Type"] == "Float":
                    info_dict[entry[0]] = list(map(lambda x: float(x), entry[1].split(",")))
                else:
                    info_dict[entry[0]] = entry[1].split(",")
        samples_list = []

        for sample_string in line_list[9:]:
            sample_dict = OrderedDict()
            if sample_string == "./.":
                sample_dict["GT"] = ["./."]
            else:
                for key, value_list in zip(line_list[8].split(":"), sample_string.split(":")):
                    if metadata["FORMAT"][key]["Type"] == "Integer":
                        sample_dict[key] = list(map(lambda x: int(x), value_list.split(",")))
                    elif metadata["FORMAT"][key]["Type"] == "Float":
                        sample_dict[key] = list(map(lambda x: float(x), value_list.split(",")))
                    else:
                        sample_dict[key] = value_list.split(",")

            samples_list.append(sample_dict)
        return line_list[0], RecordVCF(position, line_list[2], line_list[3],
                                      alt_list, quality, filter_list,
                                      info_dict, samples_list, flags=flag_set)

    @staticmethod
    def _split_by_equal_sign(string):
        try:
            index = string.index("=")
        except ValueError:
            return string, None
        return string[:index], string[index+1:]

    def _split_by_comma_sign(self, string):
        return self._split_by_sign(string, sign=",")

    @staticmethod
    def _split_by_sign(string, sign=","):
        # ignores sign in "
        index_list = [-1]
        i = 1
        while (i < len(string)):
            if string[i] == "\"":
                i += 1
                while string[i] != "\"":
                    i += 1
            if string[i] == sign:
                index_list.append(i)
            i += 1
        index_list.append(len(string))
        return [string[index_list[j] + 1: index_list[j + 1]] for j in range(0, len(index_list) - 1)]

    def filter(self, expression):
        """

        :param expression: filtering expression
        :return:
        """
        # expression should be a function with one argument - record entry
        filtered_records, filtered_out_records = self.filter_records(expression)
        return CollectionVCF(metadata=self.metadata, records_dict=filtered_records,
                             header=self.header, samples=self.samples, from_file=False),\
               CollectionVCF(metadata=self.metadata, records_dict=filtered_out_records,
                             header=self.header, samples=self.samples, from_file=False)

    def filter_by_zygoty(self):
        """

        :return:
        """
        # splits sites by zygoty, site counts as heterozygote if even in one sample it it hetorozygote
        def filter_expression(record):
            for sample_dict in record.samples_list:
                zyg = sample_dict["GT"][0].split("/")
                if zyg[0] != zyg[1]:
                    return False
            return True

        return self.filter(filter_expression)

    def record_coordinates(self, black_list=[], white_list=[]):
        """

        :param black_list:
        :param white_list:
        :return:
        """
        coord_dict = {}
        for scaffold in self.records:
            for record in self.records[scaffold]:
                if black_list and (scaffold in black_list):
                    continue
                if white_list and (scaffold not in white_list):
                    continue
                if scaffold not in coord_dict:
                    coord_dict[scaffold] = [record.pos]
                else:
                    coord_dict[scaffold].append(record.pos)
        for scaffold in coord_dict:
            coord_dict[scaffold] = np.array(coord_dict[scaffold])
        return coord_dict

    def get_positions(self):
        """

        :return:
        """
        """
        Extracts position coordinates of records in collection
        dict of coordinates, keys are scaffolds
        """
        positions_dict = OrderedDict({})
        for scaffold in self.records:
            positions_dict[scaffold] = np.array([[record.pos] for record in self.records[scaffold]])
        return positions_dict

    def check_by_ref_and_alt(self, ref_alt_list, flag):
        """

        :param ref_alt_list:
        :param flag:
        :return:
        """
        for record in self:
            record.check_ref_alt_list(ref_alt_list, flag)

    def filter_by_ref_and_alt(self, ref_alt_list):
        """

        :param ref_alt_list:
        :return:
        """
        # structure of ref_alt_list:  [[ref1,[alt1.1, alt1.M1]], ..., [refN,[altN.1, ..., altN.MN]]]
        return self.filter(lambda record: (record.ref, record.alt_list) in ref_alt_list)

    def set_filter(self, expression, filter_name):
        """
        Sets filter_name in FILTER field if expression returns True
        :param expression:
        :param filter_name:
        :return:
        """
        for record in self:
            if expression:
                if "PASS" in record.filter_list or "." in record.filter_list:
                    record.filter_list = [filter_name]
                else:
                    record.filter_list.append(filter_name)

    def split_by_scaffolds(self):
        """

        :return:
        """
        return [CollectionVCF(metadata=self.metadata, records_dict={scaffold: self.records[scaffold]},
                              header=self.header, samples=self.samples, from_file=False)
                for scaffold in self.records]

    @staticmethod
    def _reference(record):
        nucleotides = ["A", "C", "G", "T"]
        if record.ref in nucleotides:
            return record.ref
        return "INDEL"

    def rainfall_plot(self, plot_name, base_colors=[], single_fig=True, dpi=300, figsize=(40, 40), facecolor="#D6D6D6",
                      ref_genome=None, masked_regions=None, min_gap_length=10, draw_gaps=False, suptitle=None,
                      gaps_color="#777777", masked_regions_color="#aaaaaa", logbase=2,
                      extension_list=["svg", "eps", "pdf", "png", "jpg"]):
        """

        :param plot_name:
        :param base_colors:
        :param single_fig:
        :param dpi:
        :param figsize:
        :param facecolor:
        :param ref_genome:
        :param masked_regions:
        :param min_gap_length:
        :param draw_gaps:
        :param suptitle:
        :param gaps_color:
        :param masked_regions_color:
        :param logbase:
        :param extension_list:
        :return:
        """
        # TODO: add multithreading drawing if possible and multipicture drawing
        print("Drawing rainfall plot...")
        plot_dir = "rainfall_plot"
        reference_colors = {"A": "#FBFD2B",    # yellow
                            "C": "#FF000F",     # red
                            "G": "#000FFF",     # blue
                            "T": "#4ED53F",     # green
                            "INDEL": "#000000"  # black
                            }
        if base_colors:
            reference_colors = base_colors

        num_of_regions = self.number_of_scaffolds
        positions_dict = self.get_positions()
        distances_dict = {}
        region_reference_dict = {}
        os.system("mkdir -p %s" % plot_dir)
        if single_fig:
            fig = plt.figure(1, dpi=dpi, figsize=figsize, facecolor=facecolor)
            fig.suptitle(suptitle if suptitle else "Rainfall plot", fontsize=40, fontweight='bold', y=0.94)
            sub_plot_dict = OrderedDict({})
        index = 1

        for region in self.records:
            # np.ediff1d return differences between consecutive elements in array, then 0 is added to the beginning
            distances_dict[region] = np.insert(np.ediff1d(positions_dict[region]), 0, 0)
            region_reference_dict[region] = OrderedDict({"A": [[], []],
                                                         "C": [[], []],
                                                         "G": [[], []],
                                                         "T": [[], []],
                                                         "INDEL": [[], []]})
            for i in range(0, len(self.records[region])):
                region_reference_dict[region][self._reference(self.records[region][i])][0].append(positions_dict[region][i])
                region_reference_dict[region][self._reference(self.records[region][i])][1].append(distances_dict[region][i])
            if single_fig:
                if not sub_plot_dict:
                    sub_plot_dict[region] = plt.subplot(num_of_regions, 1, index, axisbg=facecolor)
                else:
                    keys = list(sub_plot_dict.keys())
                    sub_plot_dict[region] = plt.subplot(num_of_regions, 1, index,
                                                        sharex=sub_plot_dict[keys[0]],
                                                        sharey=sub_plot_dict[keys[0]],
                                                        axisbg=facecolor)

                index += 1
                if draw_gaps:
                    if ref_genome:
                        for gap in ref_genome.gaps_dict[region]:
                            plt.gca().add_patch(plt.Rectangle((gap.location.start, 1),
                                                              gap.location.end - gap.location.start,
                                                              1024*32, facecolor=gaps_color, edgecolor='none'))
                # masked regions should be SeqRecord dict
                if masked_regions:
                    for feature in masked_regions[region].features:
                        plt.gca().add_patch(plt.Rectangle((int(feature.location.start)+1, 1),
                                                           feature.location.end - feature.location.start,
                                                           1024*32, facecolor=masked_regions_color, edgecolor='none'))

                for reference in region_reference_dict[region]:
                    plt.plot(region_reference_dict[region][reference][0],
                             region_reference_dict[region][reference][1],
                             color=reference_colors[reference],
                             marker='.', linestyle='None', label=reference)

                plt.text(-0.08, 0.5, region, rotation=0, fontweight="bold", transform=sub_plot_dict[region].transAxes,
                         fontsize=30,
                         horizontalalignment='center',
                         verticalalignment='center')
                plt.ylabel("Distanse")
                plt.axhline(y=100, color="#000000")
                plt.axhline(y=1000, color="#000000")
                plt.axhline(y=500, color="purple")
                plt.axhline(y=10, color="#000000")

        if single_fig:
            for region in sub_plot_dict:
                sub_plot_dict[region].set_yscale('log', basey=logbase)
            for extension in extension_list:
                plt.savefig("%s/%s_log_scale.%s" % (plot_dir, plot_name, extension))
            plt.close()

    def hierarchical_clustering(self, method='average', dendrogramm_max_y=2000,
                                sample_name=None, save=False, clustering_dir="clustering",
                                dendrogramm_color_threshold=1000,
                                draw_dendrogramm=True,
                                write_inconsistent=True,
                                write_correlation=True):
        """
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage
        allowed methods(used to calculate distance between clusters):
        'complete'    -   Farthest Point Algorithm
        'single'      -   Nearest Point Algorithm
        'average'     -   UPGMA algorithm, distance between clusters is calculated as average from pairwise
                           distances between elements of clusters
        'weighted     -   WPGMA algorithm
        'centroid'    -   UPGMC algorithm
        'median'      -   WPGMC algorithm
        'ward'        -   incremental algorithm
        :param method:
        :param dendrogramm_max_y:
        :param sample_name:
        :param save:
        :param clustering_dir:
        :param dendrogramm_color_threshold:
        :param draw_dendrogramm:
        :param write_inconsistent:
        :param write_correlation:
        :return:
        """
        positions_dict = self.get_positions()
        correlation_dict = OrderedDict({})
        linkage_dict = OrderedDict({})
        inconsistent_dict = OrderedDict({})
        clusters_dict = OrderedDict({})
        if draw_dendrogramm or write_correlation or write_inconsistent:
            os.system("mkdir -p %s" % clustering_dir)
        for region in positions_dict:
            #print positions_dict[region]
            distance_matrix = pdist(positions_dict[region])
            #print(distance_matrix)
            linkage_dict[region] = linkage(distance_matrix, method=method)
            if draw_dendrogramm:
                plt.figure(1, dpi=150, figsize=(50, 20))
                dendrogram(linkage_dict[region],
                           color_threshold=dendrogramm_color_threshold,
                           leaf_font_size=4,
                           distance_sort=True)
                plt.ylim(ymax=dendrogramm_max_y)
                plt.axhline(y=500, color="purple")
                plt.axhline(y=1000, color="black")
                plt.savefig("%s/clustering_%s.svg" % (clustering_dir, region))
                plt.close()

            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.cophenet.html#scipy.cluster.hierarchy.cophenet
            # calculates cophenetic correlation coefficient to estimate accuracy of clustering
            correlation_dict[region] = cophenet(linkage_dict[region], distance_matrix)[0]

            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.inconsistent.html#scipy.cluster.hierarchy.inconsistent
            # calculates inconsistent coeff

            inconsistent_dict[region] = inconsistent(linkage_dict[region])
            if write_inconsistent:
                np.savetxt("%s/inconsistent_coefficient_%s.t" % (clustering_dir, region), inconsistent_dict[region])

            # clusters_dict[region] = fcluster(linkage_dict[region], 1)
            # np.savetxt("clustering/clusters_%s.t" % region, clusters_dict[region], fmt="%i")
        if write_correlation:
            sample = sample_name
            if not sample:
                sample = self.samples[0]
            with open("%s/correlation.t" % clustering_dir, "w") as cor_fd:
                cor_fd.write("sample\t%s\n" % ("\t".join(positions_dict.keys())))
                cor_fd.write("%s\t%s\n" % (sample,
                                           "\t".join([str(correlation_dict[region]) for region in positions_dict])))

        if save:
            self.linkage_dict = linkage_dict

        return linkage_dict

    def get_clusters(self,
                     extracting_method="inconsistent",
                     threshold=0.8,
                     cluster_distance='average',
                     dendrogramm_max_y=2000,
                     sample_name=None,
                     save_clustering=False,
                     clustering_dir="clustering",
                     dendrogramm_color_threshold=1000,
                     draw_dendrogramm=True,
                     return_collection=True,
                     write_inconsistent=True,
                     write_correlation=True):

        from MACE.Parsers.CCF import RecordCCF, CollectionCCF, MetadataCCF, HeaderCCF
        if self.linkage_dict:
            linkage_dict = self.linkage_dict
        else:
            linkage_dict = self.hierarchical_clustering(method=cluster_distance,
                                                        dendrogramm_max_y=dendrogramm_max_y,
                                                        sample_name=sample_name,
                                                        save=save_clustering,
                                                        clustering_dir=clustering_dir,
                                                        dendrogramm_color_threshold=dendrogramm_color_threshold,
                                                        draw_dendrogramm=draw_dendrogramm,
                                                        write_correlation=write_correlation,
                                                        write_inconsistent=write_inconsistent)
        mut_clusters_dict = OrderedDict({})
        clusters = OrderedDict()
        for region in linkage_dict:
            clusters[region] = fcluster(linkage_dict[region], threshold, criterion=extracting_method)

        if return_collection:
            record_ccf_dict = OrderedDict()
            for region in self.records:
                # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html#scipy.cluster.hierarchy.fcluster
                clusters_dict = OrderedDict({})

                for i in range(0, len(clusters[region])):
                    if clusters[region][i] not in clusters_dict:
                        clusters_dict[clusters[region][i]] = [self.records[region][i]]
                    else:
                        clusters_dict[clusters[region][i]].append(self.records[region][i])



                record_ccf_dict[region] = [RecordCCF(collection_vcf=CollectionVCF(records_dict=dict([(region, clusters_dict[cluster])]), from_file=False),
                                                     from_records=True) for cluster in clusters_dict]

            return CollectionCCF(records_dict=record_ccf_dict, metadata=MetadataCCF(self.samples,
                                                                                    vcf_metadata=self.metadata,
                                                                                    vcf_header=self.header),
                                 header=HeaderCCF("CLUSTER_ID\tCHROM\tSTART\tEND\tDESCRIPTION".split("\t")))
        else:
            return clusters

    # methods for sum of two CollectionsVCF: no check for intersections(!!!!!!!!)
    def __add__(self, other):
        new_records_dict = deepcopy(self.records)
        for scaffold in other.scaffold_list:
            if scaffold in self.scaffold_list:
                new_records_dict[scaffold] += other.records[scaffold]
            else:
                new_records_dict[scaffold] = other.records[scaffold]
        return CollectionVCF(metadata=self.metadata, records_dict=new_records_dict,
                             header=self.header, samples=self.samples, from_file=False)

    def __radd__(self, other):
        return self.__add__(other)

    def test_thresholds(self,
                        extracting_method="inconsistent",
                        threshold=None,
                        cluster_distance='average',
                        dendrogramm_max_y=2000,
                        sample_name=None,
                        save_clustering=False,
                        testing_dir="threshold_test",
                        count_singletons=True,
                        scaffold_prefix="Region"):
        # TODO: adjust parameters of figure
        # threshold is tuple(list) of three variables: min, max, number

        # extracting_method possible values
        #   inconsistent
        #   distance
        #   maxclust
        #   monocrit
        #   monocrit

        if self.linkage_dict:
            linkage_dict = self.linkage_dict
        else:
            linkage_dict = self.hierarchical_clustering(method=cluster_distance,
                                                        dendrogramm_max_y=dendrogramm_max_y,
                                                        sample_name=sample_name,
                                                        save=save_clustering,
                                                        clustering_dir=testing_dir)

        num_of_regions = len(list(linkage_dict.keys()))

        side = int(sqrt(num_of_regions))
        if side*side != num_of_regions:
            side += 1
        sub_plot_dict = OrderedDict({})
        fig = plt.figure(2, dpi=150, figsize=(30, 30))
        fig.suptitle("Relashionship between number of clusters and threshold of %s" % extracting_method, fontsize=20)

        thresholds = threshold
        if extracting_method == "inconsistent":
            if not threshold:
                thresholds = (0.5, 1.5, 21)

        index = 1
        for region in linkage_dict:
            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html#scipy.cluster.hierarchy.fcluster
            n_clusters_list = []
            n_nonsingleton_clusters = []
            n_multiclusters = []
            n_five_plus_clusters = []
            coef_threshold_list = np.linspace(*thresholds)  # best variant 0.5, 1.5, 21
            for i in coef_threshold_list:
                clusters = fcluster(linkage_dict[region], i, criterion=extracting_method)
                n_clusters_list.append(max(clusters))

                # counting clusters with 2+ and 3+ clusters
                ua, uind = np.unique(clusters, return_inverse=True)
                counted = np.bincount(uind)
                #j = 0
                nonsingleton = 0
                multicluster = 0  # 3+
                five_plus_clusters = 0 # 5+
                for k in counted:
                    if k > 1:
                        nonsingleton += 1
                    if k > 2:
                        multicluster += 1
                    if k > 4:
                        five_plus_clusters += 1
                n_nonsingleton_clusters.append(nonsingleton)
                n_multiclusters.append(multicluster)
                n_five_plus_clusters.append(five_plus_clusters)
            sub_plot_dict[region] = plt.subplot(side, side, index, axisbg="#D6D6D6")
            #ax = plt.gca()
            #ax.set_xticks(np.arange(0.5, 2.2, 0.1))

            plt.grid()
            if count_singletons:
                plt.plot(coef_threshold_list, n_clusters_list, label="all")
            plt.plot(coef_threshold_list, n_nonsingleton_clusters, "green", label="2+")
            plt.plot(coef_threshold_list, n_multiclusters, "red", label="3+")
            plt.plot(coef_threshold_list, n_five_plus_clusters, "black", label="5+")
            plt.title("%s %s" % (scaffold_prefix, region))
            plt.legend(loc='upper right')
            plt.ylabel("Number of clusters")
            plt.xlabel("Threshold")
            #plt.axvline(x=0.8, color="purple")
            #plt.axvline(x=1.1, color="purple")

            plt.ylim(ymin=0)
            index += 1
        plt.savefig("%s/clusters_%s.svg" % (testing_dir, extracting_method))
        plt.close()

    def add_info(self, metadata_line, expression, info_name, info_value=None):
        self.metadata.add_metadata(metadata_line)
        for record in self:
            if expression(record):
                value = info_value if isinstance(info_value, list) else [] if info_value is None else [info_value]
                if info_name in record.info_dict:
                    record.info_dict[info_name] += value
                else:
                    record.info_dict[info_name] = value

    def parse_snpeff_info_record(self, string):
        effect, parameters = string.split("(")
        # remove closing bracket and split
        parameters = parameters[:-1].split("|")
        return [effect] + parameters

    def extract_snpeff_info(self, output_file):

        snpeff_info_dict_keys = "EFF", "LOS", "NMD"
        record_header_list = ["Chrom", "Pos", "Ref", "Alt"]
        snpeff_header_list = ["Effect", "Effect_Impact", "Functional_Class", "Codon_Change", "Amino_Acid_Change",
                  "Amino_Acid_Length", "Gene_Name", "Transcript_BioType", "Gene_Coding", "Transcript_ID",
                  "Exon_Rank", "Genotype_Number", "ERRORS", "WARNINGS"
                  ]

        #print(output_file)
        with open(output_file, "w") as out_fd:
            header_string = "#" + "\t".join(record_header_list + snpeff_header_list) + "\n"
            out_fd.write(header_string)
            for scaffold in self.records:
                for record in self.records[scaffold]:
                    common_part = "%s\t%i\t%s\t%s" % (scaffold, record.pos, record.ref, ",".join(record.alt_list))
                    for effect in record.info_dict["EFF"]:
                        effect_parameters = self.parse_snpeff_info_record(effect)
                        num_parameters = len(effect_parameters)
                        for i in range(0, num_parameters):
                            if effect_parameters[i] == "":
                                effect_parameters[i] = "."
                        if num_parameters < 14:
                            effect_parameters += ["." for i in range(num_parameters, 14)]
                        out_fd.write(common_part + "\t" + "\t".join(effect_parameters) + "\n")

    def count_strandness(self, prefix):
        count_dict = OrderedDict({})

        for scaffold in self.records:
            count_dict[scaffold] = np.zeros((2, 4), dtype=int)
        hor_coord_dict = {"C": 0, "G": 1}
        ver_coord_dict = {"N": 0, "P": 1, "M": 2, "B": 3}

        for scaffold in self.records:
            for record in self.records[scaffold]:
                count_dict[scaffold][hor_coord_dict[record.ref]][ver_coord_dict[record.info_dict["Fstrand"][0]]] += 1

        count_dict["all"] = sum(count_dict.values())

        for chromosome in count_dict:
            with open("%s_%s.t" % (prefix, chromosome), "w") as out_fd:
                out_list = count_dict[chromosome].tolist()
                for index, name in zip(range(0, len(out_list)), ["C", "G"]):
                    out_list[index].insert(0, name)
                out_list.insert(0, [".", "N", "P", "M", "B"])
                for string_list in out_list:
                    out_fd.write("\t".join([str(x) for x in string_list]) + "\n")

        return count_dict

    def variants_start_end(self, left, right, record_dict, skip_genes_without_five_utr=False,
                           min_five_utr_len=10):

        gene_variants_positions = []
        all_variant_start_positions = []
        all_variant_end_positions = []

        for record_id in record_dict:
            for feature in record_dict[record_id].features:
                if feature.type != "gene":
                    continue
                if skip_genes_without_five_utr:
                    for sub_feature in feature.sub_features:
                        if sub_feature.type == "five_prime_UTR" and len(sub_feature) >= min_five_utr_len:
                            break
                    else:
                        continue
                #print(feature.sub_features)
                for sub_feature in feature.sub_features:
                    if sub_feature.type != "CDS":
                        continue
                    chrom = record_id
                    strand = sub_feature.strand
                    CDS_start = sub_feature.location.start + 1 if strand == +1 else sub_feature.location.end
                    CDS_end = sub_feature.location.end if strand == +1 else sub_feature.location.start + 1
                    #region_start = CDS_start - (args.left * strand)
                    #region_end = CDS_start + (args.right * strand)

                    region_start_start = CDS_start - left if strand == +1 else CDS_start - right
                    region_start_end = CDS_start + right if strand == +1 else CDS_start + left

                    region_end_start = CDS_end - left if strand == +1 else CDS_end - right
                    region_end_end = CDS_end + right if strand == +1 else CDS_end + left
                    #print("aaa")
                    start_coordinates = []
                    end_coordinates = []
                    for variant in self.records[record_id]:
                        if region_start_start <= variant.pos <= region_start_end:
                            start_coordinates.append((variant.pos - CDS_start) * strand)
                        if region_end_start <= variant.pos <= region_end_end:
                            end_coordinates.append((variant.pos - CDS_end) * strand)
                    all_variant_start_positions += start_coordinates
                    all_variant_end_positions += end_coordinates
                    #print(feature.qualifiers)
                    gene_variants_positions.append([feature.qualifiers["Name"], strand, chrom, region_start_start,
                                                    region_start_end, start_coordinates,
                                                    region_end_start, region_end_end,
                                                    end_coordinates])
        return all_variant_start_positions, all_variant_end_positions, gene_variants_positions

    def find_location(self, record_dict, key="Ftype", strand_key="Fstrand", genes_key="Genes", genes_strand_key="Gstrand",
                      feature_type_black_list=[],
                      use_synonym=False, synonym_dict=None, add_intergenic_label=True):

        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Types of features\">" % key)
        self.metadata.add_metadata("##INFO=<ID=%s,Number=1,Type=String,Description=\"Strand of features\">" % strand_key)
        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Names of genes\">" % genes_key)
        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Strands of genes\">" % genes_strand_key)
        for record in self:
            record.find_location(record_dict, key=key, strand_key=strand_key,
                                 genes_key=genes_key, genes_strand_key=genes_strand_key,
                                 feature_type_black_list=feature_type_black_list,
                                 use_synonym=use_synonym, synonym_dict=synonym_dict,
                                 add_intergenic_label=add_intergenic_label)


class ReferenceGenome(object):
    """
    ReferenceGenome class
    """
    def __init__(self, ref_gen_file, masked_regions=None, index_file="refgen.idx", filetype="fasta", mode="index_db"):
        """

        :param ref_gen_file:
        :param masked_regions:
        :param index_file:
        :param filetype:
        :param mode: mode of parsing reference genome. Allowed variants - to_dict, index, index_db
        :return:
        """
        self.ref_gen_file = ref_gen_file
        if mode == "index_db":
            self.reference_genome = SeqIO.index_db(index_file, [ref_gen_file], filetype)
        elif mode == "index":
            self.reference_genome = SeqIO.index(ref_gen_file, filetype)
        else:
            self.reference_genome = SeqIO.to_dict(SeqIO.parse(ref_gen_file, filetype))

        self.region_list = list(self.reference_genome.keys())
        self.length = 0
        self.feature_dict = OrderedDict()
        self.region_length = OrderedDict()
        for region in self.reference_genome:
            length = len(self.reference_genome[region])
            self.length += length
            self.region_length[region] = length

        self.region_index = self.rec_index()
        self.gaps_dict = OrderedDict()
        self.masked_regions = masked_regions if masked_regions else OrderedDict()

    def __len__(self):
        return self.length

    def rec_index(self):
        index_dict = OrderedDict({})
        index_dict[self.region_list[0]] = [0, self.region_length[self.region_list[0]] - 1]
        for index in range(1, self.number_of_regions()):
            index_dict[self.region_list[index]] = [index_dict[self.region_list[index-1]][1] + 1,
                                            index_dict[self.region_list[index-1]][1] + self.region_length[self.region_list[index]]]
        return index_dict

    def get_position(self, coordinate):
        # coordinate should be 0-based
        tmp_coordinate = self.region_index[self.region_list[-1]][1] + coordinate + 1 if coordinate < 0 else coordinate
        if tmp_coordinate < 0:
            raise ValueError("Coordinate %i is too small for this genome" % tmp_coordinate)
        for region in self.region_list:
            start, end = self.region_index[region]
            if start <= tmp_coordinate <= end:
                coordinate_region = region
                break
        else:
            raise ValueError("Coordinate %i is too large for this genome" % tmp_coordinate)
        shift = tmp_coordinate - start
        return coordinate_region, shift

    def number_of_regions(self):
        return len(self.region_length)

    def find_gaps(self):
        """
        Finds gaps (N) in reference genome and writes them as SeqFeatures to self.gaps_dict.
        Keys of dict are region names.
        :return: None
        """
        gap_reg_exp = re.compile("N+", re.IGNORECASE)
        for region in self.reference_genome:
            self.gaps_dict[region] = []
            gaps = gap_reg_exp.finditer(str(self.reference_genome[region].seq))  # iterator with
            for match in gaps:
                self.gaps_dict[region].append(SeqFeature(FeatureLocation(match.start(), match.end()),
                                                         type="gap", strand=None))

    def generate_snp_set(self, size, substitution_dict=None, zygoty="homo", out_vcf="synthetic.vcf"):
        multiplier = (5 - len(substitution_dict)) if substitution_dict else 1
        unique_set = np.array([], dtype=np.int64)
        print("Generating...")

        index = 1
        while len(unique_set) < size:
            print("Iteration %i..." % index)
            print("    Generating raw set %i..." % index)
            raw_set = np.random.random_integers(0, high=self.length-1, size=size*multiplier)
            print("    Generated %i..." % len(raw_set))
            print("    Removing duplicates from raw set %i..." % index)
            raw_set = np.hstack((unique_set, raw_set))
            #print(raw_set)
            unique_set = np.unique(raw_set)
            print("    After removing  duplicates left %i.." % len(unique_set))
            #print(unique_set)
            print("    Checking filtered set %i..." % index)
            mutation_count = 0
            report_count = 1
            if substitution_dict is not None:
                tmp_list = []
                references = list(substitution_dict.keys())
                for coordinate in unique_set:
                    region, interregion_pos = self.get_position(coordinate)
                    if self.reference_genome[region][interregion_pos] in references:
                        if self.masked_regions and (region in self.masked_regions):
                            for masked_region in self.masked_regions[region].features:
                                #print(self.masked_regions[region].features)
                                #print(masked_region)
                                if interregion_pos in masked_region:
                                    break
                            else:
                                tmp_list.append(coordinate)
                                mutation_count += 1
                        else:
                            tmp_list.append(coordinate)
                            mutation_count += 1
                    if mutation_count/report_count == 500:
                        print("Generated %i" % mutation_count)
                        report_count += 1
                unique_set = np.array(tmp_list)
            index += 1
        nucleotides = {"A", "C", "G", "T"}
        alt_dict = {}
        if not substitution_dict:
            for nuc in nucleotides:
                alt_dict[nuc] = list(nucleotides - {nuc})
        else:
            for ref in substitution_dict:
                if substitution_dict[ref]:
                    alt_dict[ref] = list(set(substitution_dict[ref]) - {ref})
                else:
                    alt_dict[ref] = list(nucleotides - {ref})
        print("Writing vcf...")

        with open(out_vcf, "w") as out_fd:
            out_fd.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            out_fd.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSynthetic_uniform_sample\n")
            mutation_count = 0
            for coordinate in unique_set:
                region, interregion_pos = self.get_position(coordinate)
                ref = self.reference_genome[region][interregion_pos]

                if len(substitution_dict[ref]) != 1:
                    alt = alt_dict[ref][np.random.randint(0, len(alt_dict[ref]))]
                else:
                    alt = alt_dict[ref][0]
                if zygoty == "homo":
                    zyg = "1/1"
                elif zygoty == "hetero":
                    zyg = "0/1"
                out_fd.write("%s\t%i\t.\t%s\t%s\t.\t.\t.\tGT\t%s\n" %
                             (region, interregion_pos + 1, ref, alt, zyg))
                mutation_count += 1
                if mutation_count == size:
                    break

if __name__ == "__main__":
    #workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all"
    vcf_file = "/home/mahajrod/Genetics/MACE/example_data/Lada_et_al_2015/PmCDA1_3d_SNP.vcf"
    masking_file = "/home/mahajrod/Genetics/MCTool/example_data/LAN210_v0.10m_masked_all_not_in_good_genes.gff"
    from BCBio import GFF

    masked_regions = SeqIO.to_dict(GFF.parse(masking_file))

    reference = ReferenceGenome("/home/mahajrod/Genetics/MACE/example_data/Lada_et_al_2015/LAN210_v0.10m.fasta",
                                masked_regions=masked_regions)
    print(reference.reference_genome["chrXVI"][18075])
    #reference.generate_snp_set(10000, {"C": ["T"], "G": ["A"]}, out_vcf="synthetic.vcf")
    """
    collection = CollectionVCF(from_file=True, in_file=vcf_file)
    clusters = collection.get_clusters(sample_name="EEEEEE", save_clustering=True,
                                          extracting_method="distance",
                                          threshold=1000, cluster_distance='average',
                                          dendrogramm_max_y=2000, dendrogramm_color_threshold=1000,
                                          clustering_dir="temp")
    clusters.write("tra.ccf")
    clusters.write_gff("tra.gff")
    extracted_vcf = clusters.extract_vcf()
    extracted_vcf.write("extracted.vcf")
    """

    """


    def check_location(feature, pos):
        return True if pos in feature else False

    collection.set_location_flag(masked_regions, check_location, "BR")

    for key in collection.metadata:
        print key
        print collection.metadata[key]


    reference = ReferenceGenome("/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/LAN210_v0.9m.fasta",
                                index_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/LAN210_v0.9m.idx")
    #print(reference.reference_genome)
    reference.find_gaps()
    os.chdir(workdir)
    samples_list = sorted(os.listdir("."))

    suffix = "_GATK_best_merged.vcf"

    for sample in samples_list:
        print("Handling %s" % sample)

        os.chdir(workdir)
        os.chdir(sample)
        if "alignment_LAN210_v0.9m" not in os.listdir("."):
            continue
        os.chdir("alignment_LAN210_v0.9m")

        mutations = CollectionVCF(vcf_file=sample + suffix,
                                  from_file=True)
        print("Totaly %s mutations" % len(mutations))

        mutations.test_thresholds(save_clustering=True, testing_dir="testing_theshold_inconsistent")
        #mutations.test_thresholds(extracting_method='distance', threshold=(50, 5000, 100),
        #                          testing_dir="testing_threshold")
        #mutations.rainfall_plot(sample + suffix, ref_genome=reference, draw_gaps=True)
    """