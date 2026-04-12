import os
import glob
import argparse
import textwrap

from copy import deepcopy
from pathlib import Path
from collections import OrderedDict

import pandas as pd

from RouToolPa.Parsers.STR import CollectionSTR
from RouToolPa.Parsers.BED import CollectionBED
from RouToolPa.Parsers.GFF import CollectionGFF
from RouToolPa.Parsers.BLAST import CollectionBLAST
from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.Parsers.TableIndex import CollectionTableIndex
from RouToolPa.Collections.General import SynDict, IdList


class SyntenyRoutines:

    def __init__(self):
        pass
    # ---- methods transferred from draw_synteny.py ----

    @staticmethod
    def dist_to_next_seg_from_same_chr(df):
        temp = deepcopy(df)
        temp["dist_upstream"] = df["start"].astype("Int64") - df["end"].astype("Int64").shift(1)
        temp["dist_downstream"] = df["start"].astype("Int64").shift(-1) - df["end"].astype("Int64")

        temp["dist_end_start"] = df["start"].astype("Int64") - df["end"].astype("Int64").shift(1)
        temp["dist_end_end"] = df["end"].astype("Int64") - df["end"].astype("Int64").shift(1)
        temp["embedded"] = temp["dist_end_end"] <= 0
        return temp

    def filter_isolated_short_blocks(self, bed_df, min_block_len=1000000, max_dist_between_short_blocks=3000000):
        tmp_df = bed_df.reset_index(drop=False).set_index(["scaffold", "query"]).sort_values(by=["scaffold",
                                                                                                 "start",
                                                                                                 "end",
                                                                                                 "query",
                                                                                                 "query_start",
                                                                                                 "query_end"], axis=0)
        tmp_df["target_len"] = tmp_df["end"] - tmp_df["start"]
        tmp_df["query_len"] = tmp_df["query_end"] - tmp_df["query_start"]
        #print(tmp_df.groupby(["scaffold",
        #                         "query"], sort=False, group_keys=False).apply(dist_to_next_seg_from_same_chr))
        tmp_df = tmp_df.groupby(["scaffold",
                                 "query"], sort=False, group_keys=False).apply(self.dist_to_next_seg_from_same_chr).sort_values(by=["scaffold",
                                                                                                             "start",
                                                                                                             "end",
                                                                                                             "query",
                                                                                                             "query_start",
                                                                                                             "query_end"],
                                                                                                         axis=0)  # .reset_index(level=(0,1), drop=True)
        tmp_df["embedded"].fillna(False, inplace=True)
        # recalculate distances after removal of embedded blocks
        tmp_df = tmp_df[tmp_df["embedded"] <= 0]
        tmp_df = tmp_df.groupby(["scaffold",
                                 "query"], sort=False, group_keys=False).apply(self.dist_to_next_seg_from_same_chr).sort_values(by=["scaffold",
                                                                                                             "start",
                                                                                                             "end",
                                                                                                             "query",
                                                                                                             "query_start",
                                                                                                             "query_end"],
                                                                                                         axis=0)
        tmp_df = tmp_df[~((tmp_df["target_len"] < min_block_len) & (((tmp_df["dist_upstream"] > max_dist_between_short_blocks) | tmp_df["dist_upstream"].isna()) & (
                (tmp_df["dist_downstream"] > max_dist_between_short_blocks) | (tmp_df["dist_downstream"].isna()))))]

        return tmp_df#.reset_index(level=1, drop=False)

    @staticmethod
    def merge_adjacent_blocks(df, max_dist_between_blocks=1000000):
        temp = deepcopy(df[["start", "end", "color"]]).reset_index(drop=False)
        row_list = [list(row) for row in temp.iloc[0:1, :].itertuples(index=False)]

        for row in temp.iloc[1:, :].itertuples(index=False):
            if (row[0] == row_list[-1][0]) and (row[1] == row_list[-1][1]):
                if (row[2] - row_list[-1][3]) < max_dist_between_blocks:
                    row_list[-1][3] = max(row[3], row_list[-1][3])
                else:
                    row_list.append(list(row))
            else:
                row_list.append(list(row))
        # print(row_list)
        df = pd.DataFrame.from_records(row_list, columns=["scaffold", "query", "start", "end", "color"],
                                       index=["scaffold", "query"]) # index=["scaffold", "query"]
        df["target_len"] = df["end"] - df["start"]
        return df

    @staticmethod
    def invert_coordinates_in_synteny_table(df, scaffold_list, length_df, scaffold_column, start_column, end_column, strand_column, inverted_scaffolds_label="'"):
        temp_df = deepcopy(df)
        temp_df["length_column"] = 0
        original_index_name = temp_df.index.name
        columns_list = list(temp_df.columns)
        temp_df.reset_index(drop=False, inplace=True)
        temp_df.set_index(scaffold_column, inplace=True)
        #temp_df.to_csv("tmp", sep="\t", index=True, header=True)
        for scaffold in temp_df.index.unique():
            temp_df.loc[scaffold, "length_column"] = length_df.loc[scaffold, "length"]

        temp_df.loc[temp_df.index.isin(scaffold_list), start_column], temp_df.loc[temp_df.index.isin(scaffold_list), end_column] = temp_df.loc[temp_df.index.isin(scaffold_list), "length_column"] - temp_df.loc[temp_df.index.isin(scaffold_list), end_column], \
                   temp_df.loc[temp_df.index.isin(scaffold_list), "length_column"] - temp_df.loc[temp_df.index.isin(scaffold_list), start_column]

        plus_indexes, minus_indexes = (temp_df[strand_column] == "+") & temp_df.index.isin(scaffold_list), (temp_df[strand_column] == "-") & temp_df.index.isin(scaffold_list)
        temp_df.loc[plus_indexes, strand_column], temp_df.loc[minus_indexes, strand_column] = "-", "+"
        temp_df.reset_index(drop=False, inplace=True)
        if inverted_scaffolds_label is not None:
            for scaffold in scaffold_list:
                temp_df.loc[temp_df[scaffold_column] == scaffold, scaffold_column] = scaffold + inverted_scaffolds_label
        if (original_index_name is not None) and original_index_name in temp_df.columns:
            temp_df.set_index(original_index_name, inplace=True)

        return temp_df[columns_list]  # remove added length column and restore column order

    # ---- end of methods transferred from draw_synteny.py ----
