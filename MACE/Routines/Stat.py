#!/usr/bin/env python
"""
VCF Parser Module based on pandas
"""
__author__ = 'Sergei F. Kliver'

import os
import re
import datetime

from math import sqrt
from copy import deepcopy
from functools import reduce, partial
from collections import OrderedDict
from collections.abc import Iterable

import numpy as np
import pandas as pd

#from scipy.spatial.distance import pdist # this import was moved into methods as temporary solution because of issues with scipy on rapunzel
#from scipy.cluster.hierarchy import linkage, dendrogram, inconsistent, cophenet, fcluster # this import was moved into methods as temporary solution because of issues with scipy on rapunzel

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

from RouToolPa.Collections.General import IdList, IdSet, SynDict
from RouToolPa.Parsers.GFF import CollectionGFF
from RouToolPa.Parsers.VCF import CollectionVCF
from RouToolPa.Routines import DrawingRoutines
from RouToolPa.GeneralRoutines.File import FileRoutines
import RouToolPa.Formats.VariantFormats as VariantFormats

ref_alt_variants = {"deaminases": [("C", ["T"]), ("G", ["A"])]
                    }


class StatsVCF(FileRoutines):
    """
    CollectionVCF class

    """

    def __init__(self,):
        FileRoutines.__init__(self)

    # ------------------------------------ General stats ---------------------------------------------
    @staticmethod
    def count_zygoty(collection_vcf, outfile=None):
        # suitable onl for diploid genomes
        if collection_vcf.parsing_mode in collection_vcf.parsing_modes_with_genotypes:
            zygoty_counts = OrderedDict()
            variant_number = np.shape(collection_vcf.records)[0]
            for sample in collection_vcf.samples:
                zygoty_counts[sample] = OrderedDict({
                                                     "homo": 0,
                                                     "hetero": 0,
                                                     "ref": 0,
                                                     "absent": 0,
                                                     })
                zygoty_counts[sample]["absent"] = np.sum(collection_vcf.records[sample]["GT"][0].isna() | collection_vcf.records[sample]["GT"][1].isna())
                zygoty_counts[sample]["hetero"] = np.sum(collection_vcf.records[sample]["GT"][0] != collection_vcf.records[sample]["GT"][1]) - zygoty_counts[sample]["absent"]
                zygoty_counts[sample]["ref"] = np.sum((collection_vcf.records[sample]["GT"][0] == 0) & (collection_vcf.records[sample]["GT"][1] == 0))
                zygoty_counts[sample]["homo"] = variant_number - zygoty_counts[sample]["hetero"] - zygoty_counts[sample]["absent"] - zygoty_counts[sample]["ref"]
                #collection_vcf.records.xs('GT', axis=1, level=1, drop_level=False).apply()
            zygoty_counts  = pd.DataFrame(zygoty_counts)
            if outfile:
                zygoty_counts.to_csv(outfile, sep="\t", header=True, index=True)
            return zygoty_counts
        else:
            raise ValueError("ERROR!!! Zygoty can't be counted for parsing mode used in CollectionVCF class: %s."
                             "Use 'coordinates_and_genotypes', 'genotypes' or 'complete modes'" % collection_vcf.parsing_mode)

    @staticmethod
    def count_variants(collection_vcf, outfile=None):
        if collection_vcf.parsing_mode in collection_vcf.parsing_modes_with_genotypes:
            variant_counts = OrderedDict()
            variant_number = np.shape(collection_vcf.records)[0]

            for sample in collection_vcf.samples:
                variant_counts[sample] = (collection_vcf.records[sample]["GT"].sum(axis=1) != 0).sum()

            variant_counts = pd.Series(variant_counts)
            if outfile:
                variant_counts.to_csv(outfile, sep="\t", header=True, index=True)
            return variant_counts
        else:
            raise ValueError("ERROR!!! Zygoty can't be counted for parsing mode used in CollectionVCF class: %s."
                             "Use 'coordinates_and_genotypes', 'genotypes' or 'complete' modes" % collection_vcf.parsing_mode)

    @staticmethod
    def count_singletons(collection_vcf, output_prefix=None):
        """
        Works on diploid(!!!) genomes only
        :param collection_vcf:
        :param output_prefix:
        :return:
        """
        if collection_vcf.parsing_mode in collection_vcf.parsing_modes_with_genotypes:
            singleton_counts = OrderedDict()

            double_singleton_df = []
            het_singleton_df = []
            hom_singleton_df = []
            #singleton_df = []

            for element in collection_vcf.records[collection_vcf.samples].fillna(0).itertuples(index=False, name=None):
                # ele_dict = dict(numpy.unique(a, return_counts=True))
                ele_arr = np.array(element)
                count_dict = dict(zip(*np.unique(ele_arr, return_counts=True)))
                counts = np.array([count_dict[entry] for entry in ele_arr])

                unique_mask = counts == 1
                double_singleton_df.append(unique_mask[::2] & unique_mask[1:][::2])
                het_singleton_df.append(unique_mask[::2] != unique_mask[1:][::2])
                hom_singleton_df.append((counts == 2)[::2] & (ele_arr[::2] == ele_arr[1:][::2]))  # check if allel is double and if it is homozygote
                #singleton_df.append(double_singleton_df[-1] | het_singleton_df[-1] | hom_singleton_df[-1])

            df_list = []
            singleton_counts_df = []
            for df in double_singleton_df, het_singleton_df, hom_singleton_df:#, singleton_df:
                df_list.append(pd.DataFrame.from_records(df, columns=collection_vcf.samples, index=collection_vcf.records.index))
                singleton_counts_df.append(df_list[-1].sum())

            singleton_counts_df = pd.DataFrame(singleton_counts_df, index=["double", "hetero", "homo"]).transpose()
            singleton_counts_df["all"] = singleton_counts_df.apply(sum, axis=1)
            singleton_counts_df.index.name = "sample"

            for index, typeeee in zip([0, 1, 2], ("dosi", "hesi", "hosi")):
                df_list[index].columns = pd.MultiIndex.from_arrays([collection_vcf.samples,
                                                                   [typeeee] * len(collection_vcf.samples),
                                                                   np.zeros(len(collection_vcf.samples), dtype=int)])

            merged_df = pd.concat([collection_vcf.records] + df_list, axis=1)

            """
            all_genotype_sum = (collection_vcf.records[collection_vcf.samples].xs('GT', axis=1, level=1, drop_level=False).fillna(0) != 0).sum(axis=1)
            for sample in collection_vcf.samples:
                singleton_counts[sample] = 0

                sample_genotype_sum = (collection_vcf.records[sample]["GT"].fillna(0) != 0).sum(axis=1)
                singleton_mask = (all_genotype_sum == sample_genotype_sum)
                singleton_counts[sample] = singleton_mask.sum()
                if output_prefix:
                    pass
            singleton_counts = pd.Series(singleton_counts)

            """

            if output_prefix:
                singleton_counts_df.to_csv("%s.singletons.counts" % output_prefix, sep="\t", header=True, index=True)
                merged_df["POS"] += 1
                merged_df.reset_index(level=1, drop=True, inplace=True)
                merged_df.to_csv("%s.singletons.df" % output_prefix, sep="\t", header=True, index=True)
            merged_df["POS"] -= 1

            return singleton_counts_df, merged_df
        else:
            raise ValueError("ERROR!!! Zygoty can't be counted for parsing mode used in CollectionVCF class: %s."
                             "Use 'coordinates_and_genotypes', 'genotypes' or 'complete' modes" % collection_vcf.parsing_mode)

    @staticmethod
    def extract_singletons(collection_vcf, output_prefix=None):
        if collection_vcf.parsing_mode in collection_vcf.parsing_modes_with_genotypes:
            singleton_df_dict = OrderedDict()
            variant_number = np.shape(collection_vcf.records)[0]

            all_genotype_sum = (collection_vcf.records[collection_vcf.samples].xs('GT', axis=1, level=1, drop_level=False) != 0).sum(axis=1)
            for sample in collection_vcf.samples:

                sample_genotype_sum = (collection_vcf.records[[sample]].xs('GT', axis=1, level=1, drop_level=False) != 0).sum(axis=1)

                singleton_df_dict[sample] = collection_vcf[(all_genotype_sum == sample_genotype_sum)]
            if output_prefix:
                for sample in collection_vcf.samples:
                    singleton_df_dict[sample].to_csv("%s.singletons.%s.tsv".format(output_prefix, sample), sep="\t", header=True, index=True)
            return singleton_df_dict
        else:
            raise ValueError("ERROR!!! Zygoty can't be counted for parsing mode used in CollectionVCF class: %s."
                             "Use 'coordinates_and_genotypes', 'genotypes' or 'complete' modes" % collection_vcf.parsing_mode)
    # ------------------------------------ General stats end ------------------------------------------

    # ------------------------------------ Window-based stats -----------------------------------------
    @staticmethod
    def count_window_number_in_scaffold(scaffold_length, window_size, window_step):
        if scaffold_length < window_size:
            return 0
        return int((scaffold_length - window_size)/window_step) + 1
    # ---------------------------In progress--------------------------

    def count_variants_in_windows(self, collection_vcf, window_size, window_step, reference_scaffold_lengths=None,
                                  ignore_scaffolds_shorter_than_window=True, output_prefix=None,
                                  skip_empty_windows=False, expression=None, per_sample_output=False,
                                  scaffold_black_list=None, scaffold_white_list=None,
                                  scaffold_syn_dict=None
                                  ):
        # TODO: rewrite tosimplify this function
        window_stepppp = window_size if window_step is None else window_step

        if window_stepppp > window_size:
            raise ValueError("ERROR!!! Window step(%i) can't be larger then window size(%i)" % (window_stepppp, window_size))
        elif (window_size % window_stepppp) != 0:
            raise ValueError("ERROR!!! Window size(%i) is not a multiple of window step(%i)..." % (window_size, window_stepppp))

        steps_in_window = window_size // window_stepppp
        if reference_scaffold_lengths is not None:
            tmp_len_df = reference_scaffold_lengths
        else:
            tmp_len_df = collection_vcf.metadata["contig"]

        if isinstance(tmp_len_df, pd.DataFrame):
            ref_scaf_len_df = tmp_len_df
        else:
            ref_scaf_len_df = pd.DataFrame.from_dict(tmp_len_df, orient="index")
            ref_scaf_len_df.columns = ["length"]

        def count_windows(scaffold_length):
            return self.count_window_number_in_scaffold(scaffold_length, window_size, window_stepppp)

        number_of_windows_df = ref_scaf_len_df.applymap(count_windows)
        number_of_windows_df.columns = ["WINDOW"]

        short_scaffolds_ids = number_of_windows_df[number_of_windows_df['WINDOW'] == 0]
        short_scaffolds_ids = IdSet(short_scaffolds_ids.index.unique().to_list())

        vcf_scaffolds = set(collection_vcf.scaffold_list)
        reference_scaffolds = set(ref_scaf_len_df.index.unique().to_list())

        scaffolds_absent_in_reference = IdSet(vcf_scaffolds - reference_scaffolds)

        if scaffolds_absent_in_reference:
            print(scaffolds_absent_in_reference)
            raise ValueError("ERROR!!! Some scaffolds from vcf file are absent in reference...")
        scaffolds_absent_in_vcf = IdSet(reference_scaffolds - vcf_scaffolds)

        number_of_windows_non_zero_df = number_of_windows_df[number_of_windows_df['WINDOW'] > 0]

        count_index = [[], []]
        for scaffold in number_of_windows_non_zero_df.index:
            #print("AAAa")
            #print(number_of_windows_non_zero_df.loc[scaffold])
            #print("\t{0}".format(number_of_windows_non_zero_df.loc[scaffold][0]))
            #print("\t{0}".format(number_of_windows_non_zero_df.loc[scaffold].iloc[0]))
            count_index[0] += [scaffold] * number_of_windows_non_zero_df.loc[scaffold].iloc[0]
            count_index[1] += list(np.arange(number_of_windows_non_zero_df.loc[scaffold].iloc[0]))
        count_index = pd.MultiIndex.from_arrays(count_index, names=("CHROM", "WINDOW"))

        def get_overlapping_window_indexes(step_index):
            # this function is used if windows have overlapps
            # DO NOT FORGET TO REPLACE WINDOW INDEXES EQUAL OR LARGE TO WINDOW NUMBER
            return [window_index for window_index in range(max(step_index - steps_in_window + 1, 0),
                                                           step_index + 1)]

        def convert_step_counts_to_win_counts(df, number_of_steps_per_window):

            return reduce(lambda x, y: x.add(y, fill_value=0),
                          [df.shift(periods=entry, fill_value=0) for entry in range(0, -number_of_steps_per_window, -1)])

            #if step_index < window_number else window_number)]

        if expression:
            step_index_df = collection_vcf.records[collection_vcf.records.apply(expression, axis=1)][['POS']] // window_stepppp
        else:
            step_index_df = collection_vcf.records[['POS']] // window_stepppp
        step_index_df.columns = ["WINDOWSTEP"]

        if per_sample_output:
            count_df = []

            variant_presence = collection_vcf.check_variant_presence()
            # string below is temporaly remove beforer git pull
            variant_presence.columns = collection_vcf.samples

            for sample in collection_vcf.samples:
                tmp_count_df = pd.DataFrame(0, index=count_index,
                                            columns=[sample],
                                            dtype=np.int64)
                variant_counts = step_index_df[variant_presence[sample]].reset_index(level=1).set_index(['WINDOWSTEP'], append=True).groupby(["CHROM", "WINDOWSTEP"]).count()
                tmp_count_df[tmp_count_df.index.isin(variant_counts.index)] = variant_counts
                count_df.append(tmp_count_df)
            count_df = pd.concat(count_df, axis=1)

        else:
            count_df = pd.DataFrame(0, index=count_index,
                                    columns=["All"] if len(collection_vcf.samples) > 1 else collection_vcf.samples,
                                    dtype=np.int64)
            #print("XXXXXXXX")
            #print(count_df)
            # code for staking windows: in this case window step index  is equal to window index
            variant_counts = step_index_df.reset_index(level=1).set_index(['WINDOWSTEP'],
                                                                          append=True).groupby(["CHROM",
                                                                                                "WINDOWSTEP"], group_keys=False).count()
            #print(variant_counts.index.nlevels)
            #print(variant_counts)
            count_df[count_df.index.isin(variant_counts.index)] = variant_counts
            #print("AAAAAA")
            #print(count_df)
            #bbb[bbb.index.get_level_values('CHROM').isin(number_of_windows_non_zero_df.index)]

        #print("FFFFF")
        #print(count_df)

        #print(count_df)
        #print(steps_in_window)
        if window_stepppp != window_size:
            count_df = count_df.groupby("CHROM", group_keys=False).apply(partial(convert_step_counts_to_win_counts,
                                                                         number_of_steps_per_window=steps_in_window))

            #print(count_df)
            # window_index_df = step_index_df.applymap(get_overlapping_window_indexes)
            pass

        #if ((scaffold_black_list is not None) and (not scaffold_black_list.empty)) or (scaffold_white_list is not None and (not scaffold_white_list.empty)):
        scaffold_to_keep = self.get_filtered_entry_list(count_df.index.get_level_values(level=0).unique().to_list(),
                                                        entry_black_list=scaffold_black_list,
                                                        entry_white_list=scaffold_white_list)
        count_df = count_df[count_df.index.isin(scaffold_to_keep, level=0)]

        if scaffold_syn_dict:
            count_df.rename(index=scaffold_syn_dict, inplace=True)


        # TODO: add storing of variants in uncounted tails
        #uncounted_tail_variants_number_dict = SynDict()
        #uncounted_tail_variant_number = step_size_number_df[step_size_number_df > ref_scaf_len_df].groupby(collection_vcf.records.index(level=0)).size()
        #print count_dict[collection_vcf.samples[0]][list(count_dict[collection_vcf.samples[0]].keys())[5]]
        #print "BBBBBBBBBBBBBB"
        #print count_dict[collection_vcf.samples[0]]
        if output_prefix:
            scaffolds_absent_in_reference.write("%s.scaffolds_absent_in_reference.ids" % output_prefix)
            scaffolds_absent_in_vcf.write("%s.scaffolds_absent_in_vcf.ids" % output_prefix)
            short_scaffolds_ids.write("%s.short_scaffolds.ids" % output_prefix)
            #uncounted_tail_variants_number_dict.write("%s.uncounted_tail_variant_number.tsv" % output_prefix)

            count_df.to_csv("%s.variant_counts.tsv" % output_prefix, sep='\t', header=True, index=True)

            stat_df = count_df.groupby("CHROM").agg(["mean", "median", "min", "max"])
            stat_df.to_csv("%s.variant_counts.stats" % output_prefix, sep='\t', header=True, index=True)

        #print("DDDDDDD")
        #print(count_df)
        return count_df

    @staticmethod
    def convert_variant_count_to_feature_df(count_df,  window_size, window_step, window_column="window",
                                            scaffold_column="#scaffold", value_column=None):
        # TODO: adjust this function and count_variant to merge them. Now they are separated to keep compatibility
        # function relies that there is only one track in file
        feature_df = count_df.copy()
        feature_df.index.names = ["scaffold", "window"]
        #feature_df.rename(index={"CHROM": "scaffold", "WINDOW": "window"}, inplace=True)
        feature_df.reset_index(level=1, inplace=True, drop=False)
        track_df = feature_df.copy()
        columns = feature_df.columns
        #print(feature_df)
        #print(columns)

        feature_df["start"] = feature_df["window"] * window_step
        feature_df["end"] = feature_df["start"] + window_size

        track_df["start"] = track_df["window"] * window_step
        track_df["end"] = (track_df["window"] + 1) * window_step

        return feature_df[list(columns[:-2]) + ["start", "end", "window", columns[-1] if value_column is None else value_column]], \
               track_df[list(columns[:-2]) + ["start", "end", "window", columns[-1] if value_column is None else value_column]]

    def count_feature_length_in_windows(self, collection_gff, window_size, window_step,
                                        reference_scaffold_length_df,
                                        ignore_scaffolds_shorter_than_window=True, output_prefix=None,
                                        skip_empty_windows=False, expression=None, per_sample_output=False,
                                        scaffold_black_list=(), scaffold_white_list=(),
                                        scaffold_syn_dict=None):
        window_stepppp = window_size if window_step is None else window_step

        if window_stepppp > window_size:
            raise ValueError(
                "ERROR!!! Window step(%i) can't be larger then window size(%i)" % (window_stepppp, window_size))
        elif (window_size % window_stepppp) != 0:
            raise ValueError(
                "ERROR!!! Window size(%i) is not a multiple of window step(%i)..." % (window_size, window_stepppp))

        steps_in_window = window_size / window_stepppp

        def count_windows(scaffold_length):
            return self.count_window_number_in_scaffold(scaffold_length, window_size, window_stepppp)

        number_of_windows_df = reference_scaffold_length_df.applymap(count_windows)
        number_of_windows_df.columns = ["WINDOW"]

        short_scaffolds_ids = number_of_windows_df[number_of_windows_df['WINDOW'] == 0]
        short_scaffolds_ids = IdSet(short_scaffolds_ids.index.unique().to_list())

        gff_scaffolds = set(reference_scaffold_length_df.scaffold_list)
        reference_scaffolds = set(reference_scaffold_length_df.index.unique().to_list())

        scaffolds_absent_in_reference = IdSet(gff_scaffolds - reference_scaffolds)
        if scaffolds_absent_in_reference:
            print(scaffolds_absent_in_reference)
            raise ValueError("ERROR!!! Some scaffolds from gff file are absent in reference...")
        scaffolds_absent_in_gff = IdSet(reference_scaffolds - gff_scaffolds)
        number_of_windows_non_zero_df = number_of_windows_df[number_of_windows_df['WINDOW'] > 0]

        count_index = [[], []]
        for scaffold in number_of_windows_non_zero_df.index:
            count_index[0] += [scaffold] * number_of_windows_non_zero_df.loc[scaffold][0]
            count_index[1] += list(np.arange(number_of_windows_non_zero_df.loc[scaffold][0]))
        count_index = pd.MultiIndex.from_arrays(count_index, names=("CHROM", "WINDOW"))

        count_df = pd.DataFrame(0, index=count_index, columns="counts")

        for scaffold_id in number_of_windows_non_zero_df.index: #self.region_length:
            #uncounted_tail_variants_number_dict[scaffold_id] = 0
            if scaffold_id not in collection_gff.scaffold_list:
                continue

            for (start, end) in collection_gff.records.loc[scaffold_id][["start", "end"]].itertuples(index=False):#masked_region_dict[scaffold_id]:
                max_start_step = start / window_stepppp
                min_start_step = max(max_start_step - steps_in_window + 1, 0)
                max_end_step = (end - 1) / window_stepppp
                min_end_step = max(max_end_step - steps_in_window + 1, 0)

                if min_start_step >= number_of_windows_non_zero_df.loc[scaffold_id]: #number_of_windows:
                    break

                for i in range(min_start_step, min(max_start_step + 1, min_end_step, number_of_windows_non_zero_df.loc[scaffold_id])):
                    count_df.loc[(scaffold_id, i), "counts"] += (i * window_step) + window_size - start

                for i in range(min_end_step, min(max_start_step + 1,  number_of_windows_non_zero_df.loc[scaffold_id])):
                    count_df.loc[(scaffold_id, i), "counts"] += end - start

                for i in range(max_start_step + 1, min(min_end_step,  number_of_windows_non_zero_df.loc[scaffold_id])):
                    count_df.loc[(scaffold_id, i), "counts"] += window_size

                for i in range(max(max_start_step + 1, min_end_step), min(max_end_step + 1, number_of_windows_non_zero_df.loc[scaffold_id])):
                    count_df.loc[(scaffold_id, i), "counts"] += end - (i * window_stepppp)

        if output_prefix:
            count_df.to_csv("%s.gapped_and_masked_site_counts.tsv" % output_prefix, index=True, sep="\t")

        return count_df

    # ------------------------- Distance based stats ------------------------
    @staticmethod
    def get_distances(collection_vcf):

        distance_df = deepcopy(collection_vcf.records[["POS"]])

        distance_df["distance"] = distance_df["POS"].diff()
        for scaffold in distance_df.index.get_index_level_values(level=0):
            distance_df.loc[scaffold]["distance"][0] = 0

        distance_df["distance"].astype("Int64", copy=False)

        return distance_df

    @staticmethod
    def get_linkage_for_hierarchical_clustering(vcf_df, method='average', output_prefix=None):
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
        """
        from scipy.spatial.distance import pdist
        from scipy.cluster.hierarchy import linkage, dendrogram, inconsistent, cophenet, fcluster
        per_scaffold_counts = vcf_df.groupby(level=0).count()

        vcf_df_filtered = vcf_df[["POS"]][vcf_df.index.isin(per_scaffold_counts[per_scaffold_counts["POS"] > 1].index,
                                                            level=0)]

        print("%s\tCalculating distances..." % str(datetime.datetime.now()))
        linkage_df = pd.DataFrame({"distance": vcf_df_filtered.groupby(level=0).apply(pdist)})

        print("%s\tCalculating linkage..." % str(datetime.datetime.now()))
        linkage_df["linkage"] = linkage_df["distance"].agg(linkage, method=method)

        print("%s\tCalculating inconsistency..." % str(datetime.datetime.now()))
        linkage_df["inconsistent"] = linkage_df["linkage"].agg(inconsistent)
        # https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.cluster.hierarchy.cophenet.html

        print("%s\tCalculating cophenet coefficient..." % str(datetime.datetime.now()))
        linkage_df["cophenet"] = linkage_df.apply(lambda r: cophenet(r["linkage"], r["distance"])[0], axis=1)

        if output_prefix:
            linkage_df.to_csv("%s.linkage_df.tab" % output_prefix, sep="\t", index_label="scaffold")
        
        return linkage_df

    @staticmethod
    def get_clusters(linkage_df, extracting_method="inconsistent", threshold=0.8, depth=2):
        from scipy.cluster.hierarchy import fcluster
        clusters = linkage_df["linkage"].agg(fcluster, t=threshold, criterion=extracting_method, depth=depth)

        cluster_df = []
        index = []
        for scaf in clusters.index:
            cluster_df.append(clusters[scaf])
            index += [scaf] * len(clusters[scaf])

        cluster_df = np.concatenate(cluster_df)

        return pd.DataFrame(cluster_df, index=index, columns=["clusters"])

    def test_clustering_thresholds_from_linkage(self, linkage_df,
                                                extracting_method="inconsistent",
                                                threshold_tuple=None,
                                                min_threshold=None,
                                                max_threshold=None,
                                                threshold_number=None,
                                                threshold_step=None,
                                                output_prefix=None,
                                                inconsistency_depth=2):
        # threshold is tuple(list) of three variables: min, max, number

        # extracting_method possible values
        #   inconsistent
        #   distance
        #   maxclust
        #   monocrit
        #   monocrit

        if threshold_tuple:
            threshold_list = threshold_tuple
        elif min_threshold and max_threshold:
            if threshold_number:
                threshold_list = np.linspace(min_threshold, max_threshold, threshold_number)  # best variant 0.5, 1.5, 21
            elif threshold_step:
                threshold_list = list(np.arange(min_threshold, max_threshold, threshold_step))
                threshold_list.append(max_threshold)
            else:
                raise ValueError("ERROR!!! Neither threshold step nor threshold number was set!")
        else:
            raise ValueError("ERROR!!! Neither threshold tuple nor parameters for calculation of it were set!")

        cluster_df = self.get_clusters(linkage_df, extracting_method=extracting_method,
                                       threshold=threshold_list[0])
        cluster_df.columns = [threshold_list[0]]

        for threshold in threshold_list[1:]:
            cluster_df[threshold] = self.get_clusters(linkage_df,
                                                      extracting_method=extracting_method,
                                                      threshold=threshold)["clusters"]
        cluster_number_df = cluster_df.groupby(level=0).nunique()

        if output_prefix:
            cluster_df.to_csv("%s.cluster" % output_prefix, sep="\t", index_label="scaffold")
            cluster_number_df.to_csv("%s.cluster.counts" % output_prefix, sep="\t", index_label="scaffold")
        return cluster_df

    @staticmethod
    def test_clustering_thresholds(vcf_df, method='average', output_prefix=None,
                                   extracting_method="inconsistent", threshold_tuple=None,
                                   depth=2,
                                   min_threshold=None,
                                   max_threshold=None,
                                   threshold_number=None,
                                   threshold_step=None):
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
        """
        from scipy.spatial.distance import pdist
        from scipy.cluster.hierarchy import linkage, dendrogram, inconsistent, cophenet, fcluster
        if threshold_tuple:
            threshold_list = threshold_tuple
        elif min_threshold and max_threshold:
            if threshold_number:
                threshold_list = np.linspace(min_threshold, max_threshold, threshold_number)  # best variant 0.5, 1.5, 21
            elif threshold_step:
                threshold_list = list(np.arange(min_threshold, max_threshold, threshold_step))
                threshold_list.append(max_threshold)
            else:
                raise ValueError("ERROR!!! Neither threshold step nor threshold number was set!")
        else:
            raise ValueError("ERROR!!! Neither threshold tuple nor parameters for calculation of it were set!")
        per_scaffold_counts = vcf_df.groupby(level=0).count()

        vcf_df_filtered = vcf_df[["POS"]][vcf_df.index.isin(per_scaffold_counts[per_scaffold_counts["POS"] > 1].index,
                                                            level=0)]
        cluster_df = pd.DataFrame(index=vcf_df_filtered.index, columns=threshold_list)
        cophenet_df = pd.DataFrame(index=vcf_df_filtered.index.get_level_values(level=0).unique(), columns=["cophenet"])

        for scaffold in cophenet_df.index:
            print("%s\tHandling %s..." % (str(datetime.datetime.now()), scaffold))
            print("%s\t\t%i variants" % (str(datetime.datetime.now()), len(vcf_df_filtered.loc[[scaffold]])))
            print("%s\t\tCalculating distances..." % str(datetime.datetime.now()))
            scaffold_distance = pdist(vcf_df_filtered.loc[[scaffold]])
            print("%s\t\tCalculating linkage..." % str(datetime.datetime.now()))
            scaffold_linkage = linkage(scaffold_distance, method=method)
            print("%s\t\tCalculating cophenet coefficient..." % str(datetime.datetime.now()))
            cophenet_df.loc[scaffold, "cophenet"] = cophenet(scaffold_linkage, scaffold_distance)[0]
            print("%s\t\tCalculating clusters..." % str(datetime.datetime.now()))
            for threshold in threshold_list:
                #print len(cluster_df.loc[scaffold, threshold])
                #print len(fcluster(scaffold_linkage, t=threshold, criterion=extracting_method))
                cluster_df[threshold].loc[scaffold] = fcluster(scaffold_linkage, t=threshold,
                                                               criterion=extracting_method,
                                                               depth=depth)
        cluster_number_df = cluster_df.groupby(level=0).nunique()
        if output_prefix:
            cluster_df.to_csv("%s.cluster" % output_prefix, sep="\t", index_label=True)
            cluster_number_df.to_csv("%s.cluster.counts" % output_prefix, sep="\t", index_label="scaffold")
            cophenet_df.to_csv("%s.cophenet" % output_prefix, sep="\t", index_label="scaffold")
        #print cluster_df
        return cluster_df
    # ----------------------- Distance based stats end ----------------------

    # ----------------------------Not rewritten yet--------------------------

    def rainfall_plot(self, plot_name, dpi=300, figsize=(20, 20), facecolor="#D6D6D6",
                      ref_genome=None, min_masking_length=10, suptitle=None,
                      masking_color="#777777", logbase=2,
                      extension_list=("pdf", "png"),
                      scaffold_black_list=None, scaffold_white_list=None,
                      scaffold_ordered_list=None, sort_scaffolds=False,
                      color_expression=None,
                      default_point_color='blue',
                      dot_size=None,
                      label_fontsize=None, draw_masking=False):
        """

        :param plot_name:
        :param base_colors:
        :param single_fig:
        :param dpi:
        :param figsize:
        :param facecolor:
        :param ref_genome:
        :param masked_scaffolds:
        :param min_gap_length:
        :param draw_gaps:
        :param suptitle:
        :param gaps_color:
        :param masked_scaffolds_color:
        :param logbase:
        :param extension_list:
        :param scaffold_black_list:
        :param scaffold_white_list=:
        :param scaffold_order_list=None
        :return:

        """
        # TODO: add multithreading drawing if possible and multipicture drawing
        print("Drawing rainfall plot...")
        plot_dir = "rainfall_plot"

        os.system("mkdir -p %s" % plot_dir)

        fig = plt.figure(1, dpi=dpi, figsize=figsize ) #, facecolor=facecolor)
        fig.suptitle(suptitle if suptitle else "Rainfall plot", fontweight='bold', y=0.94, fontsize=label_fontsize) #
        sub_plot_dict = OrderedDict({})
        index = 1

        final_scaffold_list = DrawingRoutines.get_filtered_scaffold_list(self.scaffold_list,
                                                                         scaffold_black_list=scaffold_black_list,
                                                                         sort_scaffolds=sort_scaffolds,
                                                                         scaffold_ordered_list=scaffold_ordered_list,
                                                                         scaffold_white_list=scaffold_white_list,
                                                                         sample_level=False)
        num_of_scaffolds = len(final_scaffold_list)
        distances_dict = OrderedDict()
        height = 0

        if (ref_genome is not None) and draw_masking:
            masking_df = ref_genome.get_merged_gaps_and_masking()
            if min_masking_length > 1:
                masking_df.remove_small_records(min_masking_length)

        for scaffold in final_scaffold_list: # self.records
            print("Handling scaffold: %s ..." % scaffold)
            distances_dict[scaffold] = self.records.loc[scaffold, "POS"].diff()
            height = max(np.max(distances_dict[scaffold]), height)
            # pandas DataFrame diff methods return differences between consecutive elements in array,
            # and first distance is NaN always, so it is replaced by 0
            distances_dict[scaffold][0] = 0
            distances_dict[scaffold].name = 'DIST'
            if color_expression:
                colors = self.records.loc[scaffold].apply(color_expression, axis=1)
                colors.name = 'COLOR'
                distances_dict[scaffold] = pd.concat([self.records.loc[scaffold, "POS"],
                                                      distances_dict[scaffold],
                                                      colors],
                                                     axis=1)
                distances_dict[scaffold] = distances_dict[scaffold].set_index('COLOR')
                color_list = colors.index.values.unique().to_list()
            else:
                distances_dict[scaffold] = pd.concat([self.records.loc[scaffold, "POS"],
                                                      distances_dict[scaffold]],
                                                     axis=1)

        length = np.max(ref_genome.seq_lengths['length']) if ref_genome is not None else np.max(self.records["POS"])

        length *= 1.1
        if length // (10 ** 9) > 2:
            def tick_formater(x, pos):
                return '%1.1f Gbp' % (x*1e-9)
        elif length // (10 ** 6) > 200:
            def tick_formater(x, pos):
                return '%.0f Mbp' % (x*1e-6)
        elif length // (10 ** 6) > 2:
            def tick_formater(x, pos):
                return '%.1f Mbp' % (x*1e-6)

        formatter = FuncFormatter(tick_formater)

        for scaffold in final_scaffold_list:
            if not sub_plot_dict:
                sub_plot_dict[scaffold] = plt.subplot(num_of_scaffolds, 1, index) #, axisbg=facecolor)

            else:
                sub_plot_dict[scaffold] = plt.subplot(num_of_scaffolds, 1, index,
                                                      sharey=sub_plot_dict[final_scaffold_list[0]])
                                                      #sharex=sub_plot_dict[keys[0]],
                                                      #)
                                                      #facecolor=facecolor)
            sub_plot_dict[scaffold].xaxis.set_major_formatter(formatter)
            index += 1

            if ref_genome is not None:
                print("\tScaffold length:%i" % ref_genome.seq_lengths.loc[scaffold])
                plt.gca().add_patch(plt.Rectangle((1, 0),
                                                  ref_genome.seq_lengths.loc[scaffold],
                                                  height, facecolor=facecolor, edgecolor='none', alpha=0.5))
                if draw_masking:
                    for masked_region in masking_df.records.loc[scaffold].itertuples(index=False):
                        plt.gca().add_patch(plt.Rectangle((masked_region[0] + 1, 1),
                                                          masked_region[1] - masked_region[0],
                                                          height, facecolor=masking_color, edgecolor='none'))

            print("Drawing scaffold: %s ..." % scaffold)

            if color_expression:
                for color in color_list:
                    plt.scatter(distances_dict[scaffold].loc[color]['POS'],
                                distances_dict[scaffold]['DIST'],
                                color=color,
                                marker='.', s=dot_size)
            else:
                #print distances_dict[scaffold]
                #print distances_dict[scaffold]['POS']
                #print distances_dict[scaffold]['DIST']
                #print "UUUUUUUU"
                plt.scatter(distances_dict[scaffold]['POS'],
                            distances_dict[scaffold]['DIST'],
                            color=default_point_color,
                            marker='.', s=dot_size)

            plt.text(-0.13, 0.5, scaffold, rotation=0, fontweight="bold", transform=sub_plot_dict[scaffold].transAxes,
                     fontsize=label_fontsize,
                     horizontalalignment='center',
                     verticalalignment='center')
            plt.ylabel("Distanse")
            #plt.axhline(y=100, color="#000000")
            #plt.axhline(y=1000, color="#000000")
            #plt.axhline(y=500, color="purple")
            #plt.axhline(y=10, color="#000000")
            sub_plot_dict[scaffold].set_yscale('log', basey=logbase)
            sub_plot_dict[scaffold].get_xaxis().set_visible(False)
            sub_plot_dict[scaffold].spines['right'].set_color('none')
            sub_plot_dict[scaffold].spines['top'].set_color('none')
            sub_plot_dict[scaffold].spines['bottom'].set_color('none')
            plt.xlim(xmin=1, xmax=length)
            #plt.ylim(ymax=height)
            #plt.tight_layout()
        #sub_plot_dict[scaffold].unshare_x_axes(sub_plot_dict[first_scaffold])
        sub_plot_dict[final_scaffold_list[-1]].get_xaxis().set_visible(True)
        sub_plot_dict[scaffold].spines['bottom'].set_color('black')
        #plt.ylim(ymax=max_distance * 1.10)
        plt.subplots_adjust(left=0.175, bottom=0.05, right=0.95, top=0.90, wspace=None, hspace=None)
        for extension in extension_list:
            plt.savefig("%s/%s_log_scale.%s" % (plot_dir, plot_name, extension))
        plt.close()

    def check_variant_presence(self, outfile=None):
        if self.parsing_mode in self.parsing_modes_with_genotypes:

            variant_presence = pd.concat([((self.records[sample]["GT"][0].notna()) & (self.records[sample]["GT"][0] != 0)) | ((self.records[sample]["GT"][1].notna()) & (self.records[sample]["GT"][1] != 0)) for sample in self.samples], axis=1)
            variant_presence.columns = self.samples
            if outfile:
                variant_presence.to_csv(outfile, sep="\t", header=True, index=True)
            return variant_presence
        else:
            raise ValueError("ERROR!!! Variant presence can't be counted for this parsing mode: %s."
                             "Use 'coordinates_and_genotypes', 'genotypes' or 'complete modes'" % self.parsing_mode)

    def get_uniq_variants(self, output_prefix):
        variant_presence = self.check_variant_presence(outfile="%s.variant_presence" % output_prefix)
        return variant_presence[variant_presence.apply(lambda s: True if np.sum(s) == 1 else False, axis=1)]

    def count_uniq_variants(self, output_prefix, extension_list=("png",), figsize=(5, 5), dpi=200,
                            title="Unique variants"):
        if self.parsing_mode in self.parsing_modes_with_genotypes:

            variant_presence = pd.concat([((self.records[sample]["GT"][0].notna()) & (self.records[sample]["GT"][0] != 0)) | ((self.records[sample]["GT"][1].notna()) & (self.records[sample]["GT"][1] != 0)) for sample in self.samples], axis=1)
            variant_presence.columns = self.samples
            uniq_variants = variant_presence[variant_presence.apply(lambda s: True if np.sum(s) == 1 else False, axis=1)]
            uniq_variant_counts = uniq_variants.apply(np.sum)

            if output_prefix:
                #variant_presence.to_csv("%s.variant_presence" % output_prefix, sep="\t", header=True, index=True)
                #uniq_variants.to_csv("%s.uniq_variants" % output_prefix, sep="\t", header=True, index=True)
                uniq_variant_counts.to_csv("%s.uniq_variants.counts" % output_prefix, sep="\t", header=True, index=True)

            fig = plt.figure(1, figsize=figsize, dpi=dpi)

            bar_width = 0.5
            bin_coord = np.arange(len(self.samples))

            plt.bar(bin_coord, uniq_variant_counts, width=bar_width, edgecolor='white', color='blue',)

            plt.ylabel('Variants', fontweight='bold')
            plt.xlabel('Sample', fontweight='bold')
            plt.xticks(bin_coord, self.samples, rotation=45)
            plt.title(title, fontweight='bold')

            for extension in extension_list:
                plt.savefig("%s.%s" % (output_prefix, extension), bbox_inches='tight')
            plt.close()

            return uniq_variant_counts
        else:
            raise ValueError("ERROR!!! Variant presence can't be counted for this parsing mode: %s."
                             "Use 'coordinates_and_genotypes', 'genotypes' or 'complete modes'" % self.parsing_mode)

    def draw_sample_parameter_distribution(self, parameter, bin_width, output_prefix=None,
                                           extension_list=("png",), suptitle=None,
                                           xlabel=None, ylabel=None, show_median=True,
                                           show_mean=True, median_relative=False, mean_relative=False, dpi=200,
                                           subplot_size=3, xlimit=None, verbose=False, ylogbase=10):

        param = self.records.xs(parameter, axis=1, level=1, drop_level=False)
        param_mean = param.apply(np.mean)
        param_median = param.apply(np.median)
        if verbose:
            print("Median:")
            print(param_median)
            print("Mean:")
            print(param_mean)

        if median_relative:
            param = param.astype(np.float32) / param_median
            param_mean = param_mean.astype(np.float32) / param_median
            param_median = param_median.astype(np.float32) / param_median
        elif mean_relative:
            param = param.astype(np.float32) / param_mean
            param_mean = param_mean.astype(np.float32) / param_mean
            param_median = param_median.astype(np.float32) / param_mean

        param_max = param.apply(np.max)
        param_min = param.apply(np.min)

        # selection of figure size
        n = int(np.sqrt(self.sample_number))
        if n * (n + 1) >= self.sample_number:
            m = n + 1
        else:
            n += 1
            m = n
        if median_relative or mean_relative:
            if param_max > max(param_median) * 10:
                bins = np.arange(0, max(param_median) * 10, bin_width)
                bins = np.concat(bins, [max(param_max)])
            else:
                bins = np.arange(0, max(param_max), 0.1)
        else:
            print(np.max(param_median))
            print(np.max(param_median)[0])
            print(max(param_median))
            if param_max > max(param_median) * 10:
                bins = np.arange(1, max(param_median) * 10, bin_width)
                bins = np.concat(bins, [max(param_max)])
            else:
                bins = np.arange(1, max(param_max), bin_width)
        bins = np.concatenate((bins, [bins[-1] + bin_width, bins[-1] + 2 * bin_width]))

        print( "Bins:")
        print( bins)

        figure, subplot_array = plt.subplots(nrows=n, ncols=m, sharex=True, sharey=True,
                                             figsize=(m*subplot_size, n*subplot_size), dpi=dpi)
        #print subplot_array
        #print np.shape(subplot_array)
        #print n, m
        for row in range(0, n):
            for col in range(0, m):
                print( row, col)
                sample_index = row * m + col
                if ylabel and col == 0:
                    subplot_array[row][col].ylabel = ylabel
                if xlabel and row == n - 1:
                    subplot_array[row][col].xlabel = xlabel

                if sample_index >= self.sample_number:
                    continue
                sample_id = self.samples[sample_index]
                #print param[sample_id]
                # TODO: adjust function to deal not only with the first column inside parameter
                subplot_array[row][col].hist(param[sample_id][parameter][0].dropna(), bins=bins, label=sample_id)
                if show_median:
                    subplot_array[row][col].axvline(x=float(param_median[sample_id]), label="median %.2f" % param_median, color="orange")
                if show_mean:
                    subplot_array[row][col].axvline(x=float(param_mean[sample_id]), label="mean %.2f" % param_mean, color="red")
                if row == 0 and col == m - 1:
                    subplot_array[row][col].legend()
                subplot_array[row][col].set_title(sample_id)
        if suptitle:
            supt = suptitle
        elif mean_relative:
            supt = "%s distribution(Mean relative)" % parameter
        elif median_relative:
            supt = "%s distribution(Median relative)" % parameter
        else:
            supt = "%s distribution" % parameter
        plt.xlim(xmin=0)
        plt.suptitle(supt)

        if output_prefix:
            for extension in extension_list:
                plt.savefig("%s.%s" % (output_prefix, extension), bbox_inches='tight')

        xlim = xlimit if xlimit else np.max(param_median)*3
        plt.xlim(xmax=xlim, xmin=0)
        if output_prefix:
            for extension in extension_list:
                plt.savefig("%s.xlim%i.%s" % (output_prefix, xlim, extension), bbox_inches='tight')
            plt.yscale('log', basey=ylogbase)
            for extension in extension_list:
                plt.savefig("%s.xlim%i.ylog.%s" % (output_prefix, xlim, extension), bbox_inches='tight')

        plt.close()

        return param

    def get_coverage_distribution(self, output_prefix, bin_width=5, dpi=200, subplot_size=3, extension_list=("png",),
                                  verbose=False):
        if self.parsing_mode in self.parsing_modes_with_sample_coverage:
            print("Drawing coverage distribution...")
            self.draw_sample_parameter_distribution("DP", bin_width, output_prefix=output_prefix,
                                                    extension_list=extension_list,
                                                    suptitle="Coverage distribution",
                                                    xlabel="Coverage", ylabel="Variants", show_median=True,
                                                    show_mean=True, median_relative=False, mean_relative=False,
                                                    dpi=dpi, subplot_size=subplot_size,
                                                    verbose=verbose)
            print("Drawing coverage distribution relative to median...")
            self.draw_sample_parameter_distribution("DP", bin_width, output_prefix="%s.median_relative" % output_prefix,
                                                    extension_list=extension_list,
                                                    suptitle="Coverage distribution(Median relative)",
                                                    xlabel="Coverage", ylabel="Variants", show_median=True,
                                                    show_mean=True, median_relative=True, mean_relative=False,
                                                    dpi=dpi, subplot_size=subplot_size)
            print("Drawing coverage distribution relative to mean...")
            self.draw_sample_parameter_distribution("DP", bin_width, output_prefix="%s.mean_relative" % output_prefix,
                                                    extension_list=extension_list,
                                                    suptitle="Coverage distribution(Mean relative)",
                                                    xlabel="Coverage", ylabel="Variants", show_median=True,
                                                    show_mean=True, median_relative=False, mean_relative=True,
                                                    dpi=dpi, subplot_size=subplot_size)
        else:
            raise ValueError("ERROR!!! Coverage distribution can't be counted for this parsing mode: %s."
                             "Use 'pos_gt_dp' or other method parsing DP column from samples fields" % self.parsing_mode)

    def calculate_masking(self, outfile, samples=None, sample_coverage=None, min_samples=1, max_coverage=2.5, min_coverage=None):
        if self.parsing_mode in self.parsing_modes_with_sample_coverage:
            samples_to_use = samples if samples else self.samples
            coverage = self.records[samples_to_use].xs("DP", axis=1, level=1, drop_level=False)
            if sample_coverage:
                sp_coverage = pd.Series(sample_coverage, dtype=np.float32)
                sp_coverage.index = pd.MultiIndex.from_arrays([samples_to_use,
                                                               ["DP"] * len(samples_to_use),
                                                               [0] * len(samples_to_use)])
            else:
                sp_coverage = coverage.apply(np.median)
            #coverage = coverage / coverage_median
            #print sp_coverage
            #print "UUUU"
            #print coverage.apply(np.median)
            boolean_array = coverage >= (max_coverage * sp_coverage)
            if min_coverage:
                boolean_array &= coverage <= (min_coverage * sp_coverage)

            outliers = boolean_array.apply(np.sum, axis=1)
            outliers = outliers[outliers >= min_samples]
            #outliers = pd.concat([self.records[self.records.index.isin(outliers.index)]["POS"], outliers], axis=1)
            outliers = self.records[self.records.index.isin(outliers.index)]["POS"]

            print("%i variants were masked" % np.shape(outliers)[0])

            self.write_df(outliers, outfile, format="simple_bed", type="1-based")

        else:
            raise ValueError("ERROR!!! Masking can't be counted for this parsing mode: %s."
                             "Use 'pos_gt_dp' or other method parsing DP column from samples fields" % self.parsing_mode)

    #########################################################################
    #                        In progress                                    #
    #########################################################################


    #########################################################################
    # methods below were not yet rewritten for compatibility with VCFpandas #
    #########################################################################

    def no_reference_allel_and_multiallel(self, record, sample_index=None, max_allels=None):
        return record.no_reference_allel_and_multiallel(sample_index=sample_index, max_allels=max_allels)

    def filter_variants_with_reference_allel_and_multiallelic(self, sample_index=None, max_allels=None):

        def expression(record):
            return self.no_reference_allel_and_multiallel(record, sample_index=sample_index, max_allels=max_allels)

        return self.filter(expression)

    @staticmethod
    def filter_zygoty_expression(record):
            for sample_dict in record.samples_list:
                zyg = sample_dict["GT"][0].split("/")
                if zyg[0] != zyg[1]:
                    return False
            return True

    def filter_by_zygoty(self):
        """
        Splits collection based on zygoty of mutation. Mutation is counted as heterozygous even if in one sample it is hetorozygous
        :return: tuple of two CollectionVCF. First contains homozygous records, second - heterozygous
        """
        """
        def filter_expression(record):
            for sample_dict in record.samples_list:
                zyg = sample_dict["GT"][0].split("/")
                if zyg[0] != zyg[1]:
                    return False
            return True
        """
        return self.filter(self.filter_zygoty_expression)

    @staticmethod
    def filter_by_filter_presence_expression(record):
        #print record.filter_list
        for filter_entry in record.filter_list:
            #print filter_entry
            if (filter_entry != "PASS") and (filter_entry != "."):
                #print "FALSE"
                return False
        #print True
        return True

    def filter_by_filter_presence(self):
        return self.filter(self.filter_by_filter_presence_expression)

    def record_coordinates(self, black_list=[], white_list=[]):
        """
        Extracts coordinates of records in collection
        :param black_list: list of scaffolds to skip
        :param white_list: list of scaffolds to consider, other are ignored
        :return: dictionary of coordinates, keys are scaffold names, values - Numpy arrays of coordinates
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
        Extracts coordinates of records in collection
        :return: dictionary of coordinates, keys are scaffold names, values - Numpy arrays of coordinates
        """
        positions_dict = OrderedDict({})
        for scaffold in self.records:
            positions_dict[scaffold] = np.array([[record.pos] for record in self.records[scaffold]])
        return positions_dict

    def check_by_ref_and_alt(self, ref_alt_list, flag, description="No description"):
        """

        :param ref_alt_list:
        :param flag:
        :return: None
        """
        self.metadata.add_metadata("##INFO=<ID=%s,Number=0,Type=Flag,Description=\"%s\">" % (flag, description))
        for record in self:
            record.check_ref_alt_list(ref_alt_list, flag)

    def filter_by_ref_and_alt(self, ref_alt_list):
        """

        :param ref_alt_list:
        :return: None
        """
        # structure of ref_alt_list:  [[ref1,[alt1.1, alt1.M1]], ..., [refN,[altN.1, ..., altN.MN]]]
        return self.filter(lambda record: (record.ref, record.alt_list) in ref_alt_list)

    def set_filter(self, expression, filter_name):
        """
        Sets filter_name in FILTER field if expression returns True
        :param expression:
        :param filter_name:
        :return: None
        """
        for scaffold in self.records:
            for record in self.records[scaffold]:
                if expression(scaffold, record):
                    if "PASS" in record.filter_list or "." in record.filter_list:
                        record.filter_list = [filter_name]
                    else:
                        record.filter_list.append(filter_name)

    def set_filter_by_intersection_with_feature(self, annotation_dict, filter_name, mode="cross",
                                                feature_type_black_list=[]):
        """

        :param annotation_dict:
        :param filter_name:
        :param mode:
        :return:
        """
        if mode == "cross":
            def expression(scaffold, record):
                if scaffold not in annotation_dict:
                    return False
                return record.check_intersection_with_features(scaffold, annotation_dict,
                                                               feature_type_black_list=feature_type_black_list)

        elif mode == "no_cross":
            def expression(scaffold, record):
                if scaffold not in annotation_dict:
                    return True
                return not record.check_intersection_with_features(scaffold, annotation_dict,
                                                                   feature_type_black_list=feature_type_black_list)

        self.set_filter(expression, filter_name)

    def check_presence(self, chrom, position, alt_list=None):
        """
        Checks presence of variant in collection
        :param chrom:
        :param position:
        :param alt_list: optional
        :return: True if variant is present in collection(same chrom, position and optionally alts) otherwise False
        """

        if chrom not in self.scaffold_list:
            return False
        for record in self.records[chrom]:
            if record.pos > position:
                return False
            if record.pos == position:
                if alt_list:
                    if alt_list != record.alt_list:
                        return False
                return True

    def split_by_scaffolds(self):
        """

        :return:
        """
        return [CollectionVCF(metadata=self.metadata, records_dict={scaffold: self.records[scaffold]},
                              header=self.header, samples=self.samples, from_file=False)
                for scaffold in self.records]

    def get_location(self, annotation_dict, key="Loc", use_synonym=False, strand_key="strand",
                     synonym_dict=None, feature_type_black_list=[]):
        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Locations of variant\">" % key)
        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Strand\">" % strand_key)
        for scaffold in self.records:
            for record in self.records[scaffold]:
                record.get_location(scaffold, annotation_dict, key=key, use_synonym=use_synonym,
                                    synonym_dict=synonym_dict, feature_type_black_list=feature_type_black_list,
                                    strand_key=strand_key)

    @staticmethod
    def _reference(record):
        nucleotides = ["A", "C", "G", "T"]
        if record.ref in nucleotides:
            return record.ref
        return "INDEL"

    def count_heterozygous_snps(self, window_size, window_step, reference_scaffold_length_dict,
                                ignore_scaffolds_shorter_than_window=True, output_prefix=None,
                                skip_empty_windows=False, per_sample_output=False):

        def heterozygous_variant(record):
            #print record.__str__()
            #print not record.is_homozygous()
            return not record.is_homozygous()

        return self.count_variants_in_windows(window_size, window_step, reference_scaffold_length_dict,
                                              ignore_scaffolds_shorter_than_window=ignore_scaffolds_shorter_than_window,
                                              output_prefix=output_prefix,
                                              skip_empty_windows=skip_empty_windows,
                                              expression=heterozygous_variant, per_sample_output=per_sample_output)

    def draw_snps_histogram(self, window_size, window_step, output_prefix, reference_genome,
                            gaps_and_masked_positions_max_fraction=0.4,
                            expression=None, masking_gff=None, parsing_mode="parse", per_sample_output=False,
                            plot_type="concatenated",
                            xlabel="Position in genome",
                            ylabel="Number of SNPs",
                            title="SNP counts in windows",
                            suptitle="",
                            extensions=["png", ],
                            masked_or_gaped_region_mark=0,
                            figure_height_per_plot=3,
                            figure_width=12,
                            multiplier=1000):

        window_stepppp = window_size if window_step is None else window_step

        kb = multiplier / 1000
        mb = multiplier / 1000000

        normalized_ylabel = "%s per %i %s" % (ylabel, mb if mb >= 1 else kb if kb >=1 else multiplier, "Mbp" if mb >= 1 else "Kbp" if kb >= 1 else "bp")

        print( "Parsing reference and...")
        reference = ReferenceGenome(reference_genome,
                                    masked_regions=None,
                                    index_file="refgen.idx",
                                    filetype="fasta",
                                    mode=parsing_mode,
                                    black_list=[],
                                    masking_gff_list=masking_gff)
        print("Merging gaps with masking...")

        gaps_and_masked_region_window_counts = reference.count_gaped_and_masked_positions_in_windows(window_size,
                                                                                                     window_stepppp,
                                                                                                     ignore_scaffolds_shorter_than_window=True,
                                                                                                     output_prefix=output_prefix,
                                                                                                     min_gap_len=1)
        print("Counting variants in windows...")
        variant_window_counts = self.count_variants_in_windows(window_size,
                                                               window_stepppp,
                                                               reference.region_length,
                                                               ignore_scaffolds_shorter_than_window=True,
                                                               output_prefix=output_prefix,
                                                               skip_empty_windows=False,
                                                               expression=expression,
                                                               per_sample_output=per_sample_output)

        normalized_variant_window_counts = SynDict()
        filtered_normalized_variant_window_counts = SynDict()

        # normalization
        if per_sample_output:
            for sample in variant_window_counts:
                normalized_variant_window_counts[sample] = SynDict()
                filtered_normalized_variant_window_counts[sample] = SynDict()
                for scaffold_id in variant_window_counts[sample]:
                        #print sample
                        #print scaffold_id
                        #print variant_window_counts[sample][scaffold_id]
                    normalized_variant_window_counts[sample][scaffold_id] = np.divide(variant_window_counts[sample][scaffold_id].astype(float), window_stepppp - gaps_and_masked_region_window_counts[scaffold_id] + 1) * multiplier
                        #print variant_window_counts[sample][scaffold_id]
                    filtered_normalized_variant_window_counts[sample][scaffold_id] = []
        else:
            for scaffold_id in variant_window_counts:
                normalized_variant_window_counts[scaffold_id] = np.divide(variant_window_counts[scaffold_id].astype(float), window_stepppp - gaps_and_masked_region_window_counts[scaffold_id] + 1) * multiplier
                filtered_normalized_variant_window_counts[scaffold_id] = []

        # filtering
        if per_sample_output:
            for sample in variant_window_counts:
                for scaffold_id in variant_window_counts[sample]:
                    for window_index in range(0, len(variant_window_counts[sample][scaffold_id])):
                        if np.isnan(variant_window_counts[sample][scaffold_id][window_index]):
                            normalized_variant_window_counts[sample][scaffold_id][window_index] = masked_or_gaped_region_mark
                        elif float(gaps_and_masked_region_window_counts[scaffold_id][window_index])/float(window_size) > gaps_and_masked_positions_max_fraction:
                            #print variant_window_counts.keys()
                            variant_window_counts[sample][scaffold_id][window_index] = masked_or_gaped_region_mark
                            normalized_variant_window_counts[sample][scaffold_id][window_index] = masked_or_gaped_region_mark
                        else:
                            filtered_normalized_variant_window_counts[sample][scaffold_id].append(normalized_variant_window_counts[sample][scaffold_id][window_index])
        else:
            for scaffold_id in variant_window_counts:
                for window_index in range(0, len(variant_window_counts[scaffold_id])):
                    if np.isnan(variant_window_counts[scaffold_id][window_index]):
                        normalized_variant_window_counts[scaffold_id][window_index] = masked_or_gaped_region_mark
                    elif float(gaps_and_masked_region_window_counts[scaffold_id][window_index])/float(window_size) > gaps_and_masked_positions_max_fraction:
                        variant_window_counts[scaffold_id][window_index] = masked_or_gaped_region_mark #variant_window_counts[scaffold_id]
                        normalized_variant_window_counts[scaffold_id][window_index] = masked_or_gaped_region_mark
                    else:
                        filtered_normalized_variant_window_counts[scaffold_id].append(normalized_variant_window_counts[scaffold_id][window_index])
        
        if per_sample_output:
            for sample in normalized_variant_window_counts:
                normalized_variant_window_counts.write("%s.%s.normalized_variant_number.tab" % (sample, output_prefix), splited_values=True)
                filtered_normalized_variant_window_counts.write("%s.%s.filtered.normalized_variant_number.tab" % (sample, output_prefix), splited_values=True)
                
        else:
            normalized_variant_window_counts.write("%s.normalized_variant_number.tab" % output_prefix, splited_values=True)
            filtered_normalized_variant_window_counts.write("%s.filtered.normalized_variant_number.tab" % output_prefix, splited_values=True)
            
        print("Drawing...")
        if plot_type == "concatenated":
            if per_sample_output:
                data = OrderedDict()
                normalized_data = OrderedDict()
                for sample in variant_window_counts:
                    data[sample] = []
                    normalized_data[sample] = []
                    for scaffold_id in reference.region_length:
                        if scaffold_id not in variant_window_counts[sample]:
                            continue
                        len(data[sample])
                        data[sample] += list(variant_window_counts[sample][scaffold_id]) + [0, ]
                        normalized_data[sample] += list(normalized_variant_window_counts[sample][scaffold_id]) + [0, ]
                #print data
                for sample in variant_window_counts:
                    data[sample] = np.array(data[sample])
                    normalized_data[sample] = np.array(normalized_data[sample])
                    bins = np.arange(len(data[sample]))
                    #print bins
                #print data[sample]

                sample_list = list(variant_window_counts.keys())
                sample_number = len(sample_list)

                figure, subplot_list = plt.subplots(nrows=sample_number, ncols=2, sharex=True, sharey=False, figsize=(figure_width, figure_height_per_plot * sample_number))
                for row_index in range(0, sample_number):
                    if row_index > 0:
                        subplot_list[row_index][0].get_shared_x_axes().join(subplot_list[row_index][0], subplot_list[0][0])
                        subplot_list[row_index][1].get_shared_x_axes().join(subplot_list[row_index][1], subplot_list[0][1])
                    for column_index in 0, 1:

                #for subplot_index in range(0, 2 * sample_number):
                    #if subplot_index % 2 == 0:
                        if column_index == 0:
                            subplot_list[row_index][column_index].plot(bins, data[sample_list[row_index]])
                            subplot_list[row_index][column_index].set_ylabel(ylabel)
                        else:
                            subplot_list[row_index][column_index].plot(bins, normalized_data[sample_list[row_index]])
                            subplot_list[row_index][column_index].set_ylabel(normalized_ylabel)
                        subplot_list[row_index][column_index].set_xlim(xmin=0)
                        subplot_list[row_index][column_index].set_xlabel(xlabel)
                        subplot_list[row_index][column_index].set_title(self.samples[row_index])
                plt.suptitle(suptitle)
            else:
                figure, subplot_list = plt.subplots(nrows=1, ncols=2,
                                                    sharex=True, sharey=False,
                                                    figsize=(figure_width, figure_height_per_plot ))
                data = []
                normalized_data = []
                for scaffold_id in reference.region_length:
                    if scaffold_id not in variant_window_counts:
                        continue
                    data += list(variant_window_counts[scaffold_id]) + [0, ]
                    normalized_data += list(normalized_variant_window_counts[scaffold_id]) + [0, ]
                data = np.array(data)
                normalized_data = np.array(normalized_data)
                print( normalized_data)
                bins = np.arange(len(data)) #* window_step
                #print data
                #print max(data)
                #print bins
                for column_index in 0, 1:
                    if column_index == 0:
                        subplot_list[column_index].plot(bins, data)
                        subplot_list[column_index].set_ylabel(ylabel)
                    else:
                        subplot_list[column_index].plot(bins, normalized_data)
                        subplot_list[column_index].set_ylabel(normalized_ylabel)

                    subplot_list[column_index].set_xlim(xmin=0)
                    subplot_list[column_index].set_xlabel(xlabel)
                    subplot_list[column_index].set_title(title)
                plt.suptitle(suptitle)
        plt.tight_layout()
        for extension in extensions:
            plt.savefig("%s.%s" % (output_prefix, extension))

    @staticmethod
    def heterozygous_variant(record):
        #print record.__str__()
        #print not record.is_homozygous()
        return not record.is_homozygous()

    def heterozygous_sample_variant(self, record, sample_index):
        #print record.__str__()
        #print sample_index, self.samples[sample_index], not record.is_homozygous()
        return not record.is_homozygous_sample(sample_index)

    def draw_heterozygous_snps_histogram(self, window_size, window_step, output_prefix, reference_genome,
                                         gaps_and_masked_positions_max_fraction=0.4,
                                         masking_gff=None, parsing_mode="parse", per_sample_output=False,
                                         plot_type="concatenated",
                                         xlabel="Position in genome",
                                         ylabel="SNPs",
                                         title="SNP counts in windows",
                                         suptitle="",
                                         extensions=["png", ],
                                         masked_or_gaped_region_mark=0,
                                         figure_height_per_plot=3,
                                         figure_width=12,
                                         multiplier=1000):

        self.draw_snps_histogram(window_size, window_step, output_prefix, reference_genome,
                                 gaps_and_masked_positions_max_fraction=gaps_and_masked_positions_max_fraction,
                                 expression=self.heterozygous_sample_variant if per_sample_output else self.heterozygous_variant,
                                 masking_gff=masking_gff,
                                 parsing_mode=parsing_mode, per_sample_output=per_sample_output,
                                 plot_type=plot_type,
                                 xlabel=xlabel,
                                 ylabel=ylabel,
                                 title=title,
                                 suptitle=suptitle,
                                 extensions=extensions,
                                 masked_or_gaped_region_mark=masked_or_gaped_region_mark,
                                 figure_height_per_plot=figure_height_per_plot,
                                 figure_width=figure_width,
                                 multiplier=multiplier)



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

    def add_info(self, metadata_line, expression, info_name, info_value=None):
        self.metadata.add_metadata(metadata_line)
        for record in self:
            if expression(record):
                value = info_value if isinstance(info_value, list) else [] if info_value is None else [info_value]
                if info_name in record.info_dict:
                    record.info_dict[info_name] += value
                else:
                    record.info_dict[info_name] = value

    def parse_snpeff_info_record(self, string, snpeff_entry="ANN"):
        if snpeff_entry == "EFF":
            effect, parameters = string.split("(")
            # remove closing bracket and split
            parameters = parameters[:-1].split("|")
            return [effect] + parameters

        elif snpeff_entry == "ANN":
            return string.split("|")

    def extract_snpeff_info(self, output_file, snpeff_entry="ANN"):

        snpeff_info_dict_keys = "EFF", "LOS", "NMD"
        record_header_list = ["Chrom", "Pos", "Ref", "Alt", "Filter"]
        if snpeff_entry == "EFF":
            snpeff_header_list = ["Effect", "Effect_Impact", "Functional_Class", "Codon_Change", "Amino_Acid_Change",
                                  "Amino_Acid_Length", "Gene_Name", "Transcript_BioType", "Gene_Coding",
                                  "Transcript_ID", "Exon_Rank", "Genotype_Number", "ERRORS", "WARNINGS"]
        elif snpeff_entry == "ANN":
            snpeff_header_list = ["Allele", "Annotation",
                                  "Putative_impact", "Gene_Name",
                                  "Gene_ID", "Feature type",
                                  "Feature ID", "Transcript biotype",
                                  "Rank", "HGVS.c",
                                  "HGVS.p", "cDNA_position",
                                  "CDS_position", "Protein_position",
                                  "Distance_to_feature", "Errors_Warnings"
                                  ]
        else:
            raise ValueError("ERROR!!! Unknow SNPeff entry: %s. Only ANN or EFF are allowed..." % snpeff_entry)

        #print(output_file)
        with open(output_file, "w") as out_fd:
            header_string = "#" + "\t".join(record_header_list + snpeff_header_list) + "\n"
            out_fd.write(header_string)
            for scaffold in self.records:
                for record in self.records[scaffold]:
                    common_part = "%s\t%i\t%s\t%s\t%s" % (scaffold, record.pos, record.ref, ",".join(record.alt_list),
                                                          ",".join(record.filter_list))

                    if snpeff_entry not in record.info_dict:
                        continue

                    for effect in record.info_dict[snpeff_entry]:
                        effect_parameters = self.parse_snpeff_info_record(effect, snpeff_entry)
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

    def draw_info_distribution(self, info_dict_key, expression, outfile_prefix,
                               extension_list=(".svg", ".png"), bins=None,):
        scaffold_distribution = OrderedDict()
        for scaffold in self.scaffold_list:
            scaffold_distribution[scaffold] = [[], []]
            for record in self.records[scaffold]:
                #print(scaffold)
                if expression(record):
                    scaffold_distribution[scaffold][0] += record.info_dict[info_dict_key]
                else:
                    scaffold_distribution[scaffold][1] += record.info_dict[info_dict_key]

        #print(scaffold_distribution[scaffold][0])
        side = int(sqrt(self.number_of_scaffolds))
        if side*side != self.number_of_scaffolds:
            side += 1
        sub_plot_dict = OrderedDict({})
        fig = plt.figure(2, dpi=150, figsize=(15, 15))
        fig.suptitle("Distribution of %s" % info_dict_key, fontsize=20, fontweight='bold')

        index = 1
        for scaffold in self.scaffold_list:
            #print(scaffold)
            #print(scaffold_distribution[scaffold][0])
            #print(scaffold_distribution[scaffold][1])
            sub_plot_dict[scaffold] = plt.subplot(side, side, index, axisbg="#D6D6D6")
            #ax = plt.gca()
            #ax.set_xticks(np.arange(0.5, 2.2, 0.1))

            plt.grid()
            num_of_bins = bins if bins is not None else 20
            maximum = max(max(scaffold_distribution[scaffold][0]) if scaffold_distribution[scaffold][0] else 0,
                          max(scaffold_distribution[scaffold][1]) if scaffold_distribution[scaffold][1] else 0)
            if isinstance(bins, Iterable):
                if maximum > num_of_bins[-1]:
                    num_of_bins[-1] = maximum + 1
            plt.hist([scaffold_distribution[scaffold][0], scaffold_distribution[scaffold][1]], bins=num_of_bins)
            plt.xlim(xmin=0)
            plt.title("%s" % scaffold, fontweight='bold')
            #plt.legend(loc='upper right')
            plt.ylabel("Number of variants")
            plt.xlabel("%s" % info_dict_key)
            #plt.axvline(x=0.8, color="purple")
            #plt.axvline(x=1.1, color="purple")

            plt.ylim(ymin=0)
            index += 1
        plt.subplots_adjust(hspace=0.27, wspace=0.27, top=0.92, left=0.05, right=0.99, bottom=0.04)
        for extension in extension_list:
            plt.savefig("%s%s" % (outfile_prefix, extension if extension[0] == "." else ".%s" % extension))
        plt.close()

    def set_filter_for_indels_in_homopolymers(self, reference_dict, min_homopolymer_len=4,
                                              filter_name="indel_in_homopolymer"):

        def expression(scaffold, record):
            if not record.check_indel():
                False
            if scaffold not in reference_dict:
                raise ValueError("Scaffold %s is absent in reference" % scaffold)

            scaffold_length = len(reference_dict[scaffold])

            reference_pos = record.pos - 1
            left_flank = None if record.pos == 1 else reference_dict[scaffold].seq[max(0, record.pos-20):record.pos]
            right_flank = None if (record.pos + 1) == scaffold_length else reference_dict[scaffold].seq[record.pos+1: max(scaffold_length - 1,
                                                                                                                          record.pos++20)]
            ref_var_len = len(record.ref)
            indel_type_list = [None for variant in record.alt_list]

            for i in range(0, len(record.alt_list)):
                alt_variant = record.alt_list[i]
                alt_var_len = len(alt_variant)
                if len(alt_variant) < ref_var_len:
                    if record.ref[0:alt_var_len] == alt_variant:
                        indel_type_list[i] = "right"
                    elif record.ref[-alt_var_len:] == alt_variant:
                        indel_type_list[i] = "left"
                elif len(alt_variant) > ref_var_len:
                    if alt_variant[0:ref_var_len] == record.ref:
                        indel_type_list[i] = "right"
                    elif alt_variant[-ref_var_len:] == record.ref:
                        indel_type_list[i] = "left"
                else:
                    continue
                homo_len = 1
                if indel_type_list[i] == "right":
                    for letter in range(1, len(right_flank)):
                        if letter != right_flank[0]:
                            break
                        homo_len += 1
                elif indel_type_list[i] == "left":
                    for letter in range(-2, -len(right_flank)-1, -1):
                        if letter != right_flank[0]:
                            break
                        homo_len += 1
                if homo_len >= min_homopolymer_len:
                    return True
            return False
        self.set_filter(expression, filter_name)


if __name__ == "__main__":
    pass