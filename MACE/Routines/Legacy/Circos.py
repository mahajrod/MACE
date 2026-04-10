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

import numpy as np
import pandas as pd

#from scipy.spatial.distance import pdist
#from scipy.cluster.hierarchy import linkage, dendrogram, inconsistent, cophenet, fcluster

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


class Circos(FileRoutines):
    def __init__(self):
        pass
    """"
    @staticmethod
    def convert_synteny_psl_to_circos_links(collection_psl):

        links_df = collection_psl.records[[collection_psl.query_id_syn, collection_psl.query_start_syn, collection_psl.query_end_syn,
                                          collection_psl.target_id_syn, collection_psl.target_start_syn, collection_psl.target_end_syn,
                                          collection_psl.query_strand_syn]]
        links_df.loc[links_df[collection_psl.query_strand_syn] == "+-", collection_psl.target_start_syn], \
        links_df.loc[links_df[collection_psl.query_strand_syn] == "+-", collection_psl.target_end_syn] = \
                    links_df.loc[links_df[collection_psl.query_strand_syn] == "+-", collection_psl.target_end_syn], \
                    links_df.loc[links_df[collection_psl.query_strand_syn] == "+-", collection_psl.target_start_syn]

        links_df[collection_psl.query_strand_syn] = "#" + links_df[collection_psl.query_strand_syn]

        return links_df[[collection_psl.query_id_syn, collection_psl.query_start_syn, collection_psl.query_end_syn,
                         collection_psl.target_id_syn, collection_psl.target_start_syn, collection_psl.target_end_syn,
                         collection_psl.query_strand_syn]]
    """
    @staticmethod
    def swipe_end_start(row):
        return row[0], row[1] if row[2] == "+" else row[1], row[0]

    def convert_bed_synteny_track_to_links(self, collection, type="psl"):
        if type == "psl":
            links_df = collection.records[["tName", "tStart", "tEnd", "qName", "qStart", "qEnd", "strand"]],
            links_df.columns = pd.Index(["scaffold", "start", "end", "query", "query_start", "query_end", "strand"])
            links_df.set_index("scaffold", inplace=True)
        elif type == "bed":
            links_df = collection
        else:
            raise ValueError("ERROR!!! Unknown collection type")

        links_df.loc[:, ["query_start", "query_end"]] = collection.records[["query_start", "query_end", "strand"]].apply(self.swipe_end_start, axis=1, result_type='expand')

        return links_df[["start", "end", "query", "query_start", "query_end"]]

    def convert_length_df_to_karyotype(self, length_df):
        karyotype_df = length_df.reset_index(drop=False)
        karyotype_df.columns = pd.Index(["id", "end"])
        karyotype_df["-"] = "-"
        karyotype_df["chr"] = "chr"
        karyotype_df["start"] = 0
        karyotype_df["label"] = karyotype_df["id"]
        karyotype_df["color"] = karyotype_df["id"]

        return karyotype_df[["chr", "-", "id", "label", "start", "end", "color"]]
