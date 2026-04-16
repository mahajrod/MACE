import os
import glob
import argparse
import textwrap
import numpy as np
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
        # returns df with both scaffold and query as a multi level  index.
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
        tmp_df["embedded"] = tmp_df["embedded"].fillna(False)
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
    def detect_overlapping_blocks(df,
                                  query_overlapping_block_column_name,
                                  query_overlapping_fraction_column_name,
                                  query_reverse_overlapping_fraction_column_name,
                                  query_start_column_name, query_end_column_name,
                                  target_overlapping_block_column_name,
                                  target_overlapping_fraction_column_name,
                                  target_reverse_overlapping_fraction_column_name,
                                  target_start_column_name, target_end_column_name):
        output_df = deepcopy(df)
        for column_name in query_overlapping_block_column_name, target_overlapping_block_column_name:
            output_df[column_name] = pd.NA
        for column_name in query_overlapping_fraction_column_name, \
                query_reverse_overlapping_fraction_column_name, \
                target_overlapping_fraction_column_name, \
                target_reverse_overlapping_fraction_column_name:
            output_df[column_name] = np.nan

        # print(output_df)
        for row_index in range(0, len(output_df)):
            output_df["query.t_e - x_s"] = output_df[query_end_column_name].iloc[row_index] - output_df[query_start_column_name]

            output_df["query.t_e - x_s"] = output_df[query_end_column_name].iloc[row_index] - output_df[query_start_column_name]

            output_df["query.t_e - x_s"] = output_df[query_end_column_name].iloc[row_index] - output_df[query_start_column_name]

            output_df["query.t_e - x_s"] = output_df[query_end_column_name].iloc[row_index] - output_df[query_start_column_name]

        return output_df

    @staticmethod
    def detect_same_coords_blocks(df,
                                  query_same_coords_in_block_column_name,
                                  query_scaffold_id_column_name,
                                  query_start_column_name, query_end_column_name,
                                  target_same_coords_in_block_column_name,
                                  target_scaffold_id_column_name,
                                  target_start_column_name, target_end_column_name):
        output_df = deepcopy(df)
        for column_name in query_same_coords_in_block_column_name, target_same_coords_in_block_column_name:
            output_df[column_name] = pd.NA
        # print(output_df)
        for row_index in range(0, len(output_df)):
            query_same_coords_set = set(output_df[(output_df[query_scaffold_id_column_name] ==
                                                   output_df.iloc[row_index][query_scaffold_id_column_name]) & \
                                                  (output_df[query_start_column_name] == output_df.iloc[row_index][
                                                      query_start_column_name]) & \
                                                  (output_df[query_end_column_name] == output_df.iloc[row_index][
                                                      query_end_column_name])]['synteny_block_id'])
            query_same_coords_set.remove(output_df.iloc[row_index]['synteny_block_id'])

            if query_same_coords_set:
                output_df[query_same_coords_in_block_column_name].iloc[row_index] = ",".join(query_same_coords_set)

            target_same_coords_set = set(output_df[(output_df[target_scaffold_id_column_name] ==
                                                    output_df.iloc[row_index][target_scaffold_id_column_name]) & \
                                                   (output_df[target_start_column_name] == output_df.iloc[row_index][
                                                       target_start_column_name]) & \
                                                   (output_df[target_end_column_name] == output_df.iloc[row_index][
                                                       target_end_column_name])]['synteny_block_id'])
            target_same_coords_set.remove(output_df.iloc[row_index]['synteny_block_id'])

            if target_same_coords_set:
                output_df[target_same_coords_in_block_column_name].iloc[row_index] = ",".join(target_same_coords_set)

        return output_df

    @staticmethod
    def detect_nested_blocks(df,
                             nested_in_block_column_name,
                             # query_same_cooords_in_block_column_name,
                             query_nested_in_block_column_name,
                             query_scaffold_id_column_name,
                             query_start_column_name, query_end_column_name,
                             # target_same_cooords_in_block_column_name,
                             target_nested_in_block_column_name,
                             target_scaffold_id_column_name,
                             target_start_column_name, target_end_column_name):
        sorted_df = df.sort_values(by=[query_scaffold_id_column_name, query_start_column_name, query_end_column_name])
        for column_name in query_nested_in_block_column_name, target_nested_in_block_column_name:  # , query_same_cooords_in_block_column_name, target_same_cooords_in_block_column_name:
            sorted_df[column_name] = pd.NA
        # sorted_df[query_nested_in_block_column_name] = pd.NA
        # sorted_df[target_nested_in_block_column_name] = pd.NA
        query_nested_in_block_column_name_idx = sorted_df.columns.get_loc(query_nested_in_block_column_name)
        for row_index in range(0, len(sorted_df)):
            block_start = sorted_df.iloc[row_index][query_start_column_name]
            block_end = sorted_df.iloc[row_index][query_end_column_name]
            # check_df = sorted_df.iloc[:row_index]
            # print(sorted_df.iloc[row_index])
            # print(check_df)
            nested_in_block_set = set(
                sorted_df.iloc[:row_index][sorted_df.iloc[:row_index][query_end_column_name] >= block_end][
                    'synteny_block_id'])
            nested_in_block_set |= set(sorted_df.iloc[row_index + 1:][
                                           (sorted_df.iloc[row_index + 1:][query_start_column_name] == block_start) & (
                                                       sorted_df.iloc[row_index + 1:][
                                                           query_end_column_name] >= block_end)]['synteny_block_id'])

            # print(nested_in_block_set)
            if nested_in_block_set:
                # sorted_df[query_nested_in_block_column_name].iloc[row_index] = ",".join(nested_in_block_set)
                sorted_df.iloc[row_index, query_nested_in_block_column_name_idx] = ",".join(nested_in_block_set)
            # same_coords_set = set(sorted_df[(sorted_df[query_start_column_name] == block_start) & (sorted_df[query_end_column_name] == block_end)]['synteny_block_id'])
            # same_coords_set.remove(sorted_df.iloc[row_index]['synteny_block_id'])

            # if same_coords_set:
            #    sorted_df[query_same_cooords_in_block_column_name].iloc[row_index] = ",".join(same_coords_set)

        sorted_df = sorted_df.sort_values(
            by=[target_scaffold_id_column_name, target_start_column_name, target_end_column_name])
        target_nested_in_block_column_name_idx = sorted_df.columns.get_loc(target_nested_in_block_column_name)
        for row_index in range(0, len(sorted_df)):
            block_start = sorted_df.iloc[row_index][target_start_column_name]
            block_end = sorted_df.iloc[row_index][target_end_column_name]
            # check_df = sorted_df.iloc[:row_index]
            nested_in_block_set = set(
                sorted_df.iloc[:row_index][sorted_df.iloc[:row_index][target_end_column_name] >= block_end][
                    'synteny_block_id'])
            nested_in_block_set |= set(sorted_df.iloc[row_index + 1:][
                                           (sorted_df.iloc[row_index + 1:][target_start_column_name] == block_start) & (
                                                       sorted_df.iloc[row_index + 1:][
                                                           target_end_column_name] >= block_end)]['synteny_block_id'])

            if nested_in_block_set:
                sorted_df.iloc[row_index, target_nested_in_block_column_name_idx] = ",".join(nested_in_block_set)

            # same_coords_set = set(sorted_df[(sorted_df[target_start_column_name] == block_start) & (sorted_df[target_end_column_name] == block_end)]['synteny_block_id'])
            # same_coords_set.remove(sorted_df.iloc[row_index]['synteny_block_id'])
            # if same_coords_set:
            #    sorted_df[target_same_cooords_in_block_column_name].iloc[row_index] = ",".join(same_coords_set)

        def get_nested(row):
            if row.hasnans:
                return pd.NA
            else:
                return ",".join(set(row.iloc[0].split(",")) & set(row.iloc[1].split(",")))

        sorted_df[nested_in_block_column_name] = sorted_df[[query_nested_in_block_column_name, target_nested_in_block_column_name]].apply(get_nested, axis=1)
        # print(sorted_df)
        return sorted_df