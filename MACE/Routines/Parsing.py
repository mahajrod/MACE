import os
import argparse
import textwrap
from collections import OrderedDict
from copy import deepcopy

import pandas as pd

from RouToolPa.Parsers.STR import CollectionSTR
from RouToolPa.Parsers.BED import CollectionBED
from RouToolPa.Parsers.GFF import CollectionGFF
from RouToolPa.Parsers.BLAST import CollectionBLAST
from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.Parsers.TableIndex import CollectionTableIndex
from RouToolPa.Collections.General import SynDict, IdList


class ParsingRoutines:

    def __init__(self):
        """
        self.parsing_modes = {
                              "draw_features.py": {
                                    "gtf": {"class":         CollectionGFF,
                                            "format":        "gff",
                                            "parsing_mode":  "only_coordinates"},
                                    "gff": {"class":         CollectionGFF,
                                            "format":        "gff",
                                            "parsing_mode":  "only_coordinates"},
                                    "bed": {"class":         CollectionBED,
                                            "format":        "bed",
                                            "parsing_mode":  "only_coordinates"},
                                    "bedgraph": {"class":        CollectionBED,
                                                 "format":       "bed",
                                                 "parsing_mode": "only_coordinates"},
                                    "bed_colored": {"class":        CollectionBED,
                                                    "format":       "bed_colored",
                                                    "parsing_mode": "complete"},
                                    "bed_track": {"class":        CollectionBED,
                                                  "format":       "bed_track",
                                                  "parsing_mode": "all"},
                                    "bedgraph_colored": {"class": CollectionBED,
                                                         "format": "bed_track",
                                                         "parsing_mode": "all"},
                                    "bed_with_header": {"class":        CollectionBED,
                                                        "format":       "coordinates_only",
                                                        "parsing_mode": "complete"},

                                    }
                              }
        """

    @staticmethod
    def read_str_series(s):
        if s is None:
            return pd.Series()
        if os.path.exists(s):
            return pd.read_csv(s, header=None, dtype=str).squeeze("columns")
        else:
            return pd.Series(s.split(","))

    def read_mace_auxiliary_input(self,
                                  len_file=None,
                                  whitelist_file=None, max_scaffolds=50,
                                  orderlist_file=None,
                                  invertlist_file=None,
                                  queryswitchstrandlist_file=None,
                                  targetswitchstrandlist_file=None,
                                  syn_file=None, syn_file_key_column=0, syn_file_value_column=1,
                                  centromere_bed=None,
                                  highlight_bed=None,
                                  legend_file=None,
                                  vert_track_group_file=None,
                                  hor_track_group_file=None,
                                  hor_track_subgroup_file=None,
                                  masking_file=None, masking_bedgraph=None,
                                  max_masking_threshold=0.33,
                                  coverage_file=None,
                                  mean_coverage_file=None,
                                  median_coverage_file=None,
                                  min_coverage_threshold=0.33, max_coverage_threshold=2.5,
                                  median_coverage=None, mean_coverage=None):
        auxiliary_dict = OrderedDict()

        auxiliary_dict["syn_dict"] = SynDict(filename=syn_file,
                                             key_index=syn_file_key_column,
                                             value_index=syn_file_value_column)
        auxiliary_dict["max_masking_threshold"] = max_masking_threshold
        auxiliary_dict["min_coverage_threshold"] = min_coverage_threshold
        auxiliary_dict["max_coverage_threshold"] = max_coverage_threshold
        auxiliary_dict["median_coverage"] = median_coverage
        auxiliary_dict["mean_coverage"] = mean_coverage

        if centromere_bed is not None:
            if os.path.exists(centromere_bed):
                try:  # centromere_bed might be empty
                    auxiliary_dict["centromere_df"] = CollectionBED(in_file=centromere_bed, parsing_mode="coordinates_only",
                                                                    format="bed").records
                    # centromere_df.index = pd.Index(list(map(str, centromere_df.index)))
                except pd.errors.EmptyDataError:
                    auxiliary_dict["centromere_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Centromere bed file {centromere_bed} doesn't exist!")
        else:
            auxiliary_dict["centromere_df"] = None

        if len_file is not None:
            if os.path.exists(len_file):
                try:  # len_file might be empty
                    if len_file[-6:] == ".fasta":  # check if len_file is fasta
                        auxiliary_dict["len_df"] = CollectionSequence(in_file=len_file, format="fasta", parsing_mode="parse", get_stats=True).seq_lengths
                    elif len_file[-4:] == ".fai":  # check if len_file is fai
                        auxiliary_dict["len_df"] = CollectionTableIndex(in_file=len_file, format="fai", parsing_mode="length_only").records
                    else: # treat other extensions as .len file
                        auxiliary_dict["len_df"] = pd.read_csv(len_file, sep='\t', header=None, names=("scaffold", "length"),
                                                               index_col=0, dtype={"scaffold": str, "length": int})
                    #len_df.index = pd.Index(list(map(str, len_df.index)))
                except pd.errors.EmptyDataError:
                    auxiliary_dict["len_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Length file {len_file} doesn't exist!")
        else:
            auxiliary_dict["len_df"] = None

        if legend_file is not None:
            if os.path.exists(legend_file):
                try:  # legend_file might be empty
                    auxiliary_dict["legend_df"] = pd.read_csv(legend_file, header=None, index_col=0, usecols=[0, 1], sep="\t", comment=None)
                    #len_df.index = pd.Index(list(map(str, len_df.index)))
                except pd.errors.EmptyDataError:
                    auxiliary_dict["legend_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Legend file {legend_file} doesn't exist!")
        else:
            auxiliary_dict["legend_df"] = None

        if highlight_bed is not None: # TODO: reconsider the concept of the highlighting - maybe do it via special type of the track, and refactor the code
            if os.path.exists(highlight_bed):
                try:  # highlight_bed might be empty
                    auxiliary_dict["highlight_df"] = CollectionBED(in_file=highlight_bed, parsing_mode="complete",
                                                                   format="bed_colored").records #pd.read_csv(highlight_bed, header=None, index_col=0, usecols=[0, 1, 2, 3], sep="\t", comment=None)
                    #len_df.index = pd.Index(list(map(str, len_df.index)))
                except pd.errors.EmptyDataError:
                    auxiliary_dict["highlight_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Highlight bed file {highlight_bed} doesn't exist!")
        else:
            auxiliary_dict["highlight_df"] = None

        if hor_track_subgroup_file is not None:
            if os.path.exists(hor_track_subgroup_file):
                try:  # hor_track_subgroup_file might be empty
                    auxiliary_dict["hor_track_subgroup_df"] = None  # TODO: Select format (main candidate AGP) and implement parsing
                except pd.errors.EmptyDataError:
                    auxiliary_dict["hor_track_subgroup_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Horizontal track subgroup file {hor_track_subgroup_file} doesn't exist!")
        else:
            auxiliary_dict["hor_track_subgroup_df"] = None

        if vert_track_group_file is not None:
            if os.path.exists(vert_track_group_file):
                try:  # vert_track_group_file might be empty
                    auxiliary_dict["vert_track_group_df"] = None  # TODO Select format  and implement parsing
                except pd.errors.EmptyDataError:
                    auxiliary_dict["vert_track_group_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Vertical track group file {vert_track_group_file} doesn't exist!")
        else:
            auxiliary_dict["vert_track_group_df"] = None

        if hor_track_group_file is not None:
            if os.path.exists(hor_track_group_file):
                try:  # hor_track_group_file might be empty
                    auxiliary_dict["hor_track_group_df"] = None  # TODO Select format  and implement parsing
                except pd.errors.EmptyDataError:
                    auxiliary_dict["hor_track_group_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Horizontal track group file {hor_track_group_file} doesn't exist!")
        else:
            auxiliary_dict["hor_track_group_df"] = None

        if masking_file is not None:
            if os.path.exists(masking_file):
                try:  # masking_file might be empty
                    auxiliary_dict["masking_df"] = CollectionBED(in_file=masking_file, parsing_mode="coordinates_only",
                                                                 format="bed").records
                except pd.errors.EmptyDataError:
                    auxiliary_dict["masking_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Masking file {masking_file} doesn't exist!")
        else:
            auxiliary_dict["masking_df"] = None

        if masking_bedgraph is not None:
            if os.path.exists(masking_bedgraph):
                try:  # masking_file might be empty
                    auxiliary_dict["masking_bedgraph_df"] = CollectionBED(in_file=masking_bedgraph, parsing_mode="complete",
                                                                          format="bedgraph").records
                except pd.errors.EmptyDataError:
                    auxiliary_dict["masking_bedgraph_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Masking bedgraph {masking_bedgraph} doesn't exist!")
        else:
            auxiliary_dict["masking_bedgraph_df"] = None

        if coverage_file is not None:
            if os.path.exists(coverage_file):
                try:  # hor_track_group_file might be empty
                    auxiliary_dict["coverage_df"] = CollectionBED(in_file=coverage_file, parsing_mode="complete",
                                                                  format="bedgraph").records
                except pd.errors.EmptyDataError:
                    auxiliary_dict["coverage_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Coverage file {coverage_file} doesn't exist!")
        else:
            auxiliary_dict["coverage_df"] = None

        if mean_coverage_file is not None:
            if os.path.exists(mean_coverage_file):
                try:  # mean_coverage_file) might be empty
                    auxiliary_dict["mean_coverage_df"] = CollectionBED(in_file=mean_coverage_file, parsing_mode="complete",
                                                                       format="bedgraph").records
                except pd.errors.EmptyDataError:
                    auxiliary_dict["mean_coverage_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Mean coverage file {mean_coverage_file} doesn't exist!")
        else:
            auxiliary_dict["mean_coverage_df"] = None

        if median_coverage_file is not None:
            if os.path.exists(median_coverage_file):
                try:  # median_coverage_filee might be empty
                    auxiliary_dict["median_coverage_df"] = CollectionBED(in_file=median_coverage_file, parsing_mode="complete",
                                                                         format="bedgraph").records
                except pd.errors.EmptyDataError:
                    auxiliary_dict["median_coverage_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Median coverage file {median_coverage_file} doesn't exist!")
        else:
            auxiliary_dict["coverage_df"] = None

        auxiliary_dict["whitelist_series"] = self.read_str_series(whitelist_file)  # might be not only a file, but a comma-separated list too
        auxiliary_dict["orderlist_series"] = self.read_str_series(orderlist_file)  # might be not only a file, but a comma-separated list too
        auxiliary_dict["invertlist_series"] = self.read_str_series(invertlist_file)  # might be not only a file, but a comma-separated list too
        auxiliary_dict["queryswitchstrandlist_series"] = self.read_str_series(queryswitchstrandlist_file)  # might be not only a file, but a comma-separated list too
        auxiliary_dict["targetswitchstrandlist_series"] = self.read_str_series(targetswitchstrandlist_file)  # might be not only a file, but a comma-separated list too

        auxiliary_dict["max_scaffolds"] = max_scaffolds

        return auxiliary_dict

    @staticmethod
    def resolve_mace_single_genome_input(records_df, auxiliary_dict):
        if auxiliary_dict["len_df"] is None:
            raise ValueError("ERROR!!! Length dataframe was not parsed!!!")

        if auxiliary_dict["whitelist_series"].empty:
            if auxiliary_dict["orderlist_series"].empty: # Take up to auxiliary_dict["max_scaffolds"] longest scaffolds from records
                record_scaffold_len_df = auxiliary_dict["len_df"].loc[records_df.index.unique()].sort_values(by="length", ascending=False)
                auxiliary_dict["whitelist_series"] = pd.Series(record_scaffold_len_df.iloc[0:auxiliary_dict["max_scaffolds"]].index.unique())
                auxiliary_dict["orderlist_series"] = deepcopy(auxiliary_dict["whitelist_series"]).replace(auxiliary_dict["syn_dict"])
                tmp_df = records_df.loc[records_df.index.isin(auxiliary_dict["whitelist_series"])].rename(index=auxiliary_dict["syn_dict"])

            else:
                tmp_df = records_df.rename(index=auxiliary_dict["syn_dict"])
                tmp_df = tmp_df.loc[tmp_df.index.isin(auxiliary_dict["orderlist_series"])]
        else:
            tmp_df = records_df.loc[records_df.index.isin(auxiliary_dict["whitelist_series"])].rename(index=auxiliary_dict["syn_dict"])
            # remove from orderlist scaffolds which are absent in the whitelist
            if auxiliary_dict["orderlist_series"].empty:
                auxiliary_dict["orderlist_series"] = auxiliary_dict["whitelist_series"].replace(auxiliary_dict["syn_dict"])
            else:
                renamed_whitelist_series = auxiliary_dict["whitelist_series"].replace(auxiliary_dict["syn_dict"])
                # remove orderlist scaffolds that are not in whitelist
                auxiliary_dict["orderlist_series"] = auxiliary_dict["orderlist_series"].loc[auxiliary_dict["orderlist_series"].isin(renamed_whitelist_series)]
                # find whitelist scaffolds that are not in orderlist
                unordered_whitelist_series = renamed_whitelist_series[~renamed_whitelist_series.isin(auxiliary_dict["orderlist_series"])]
                # add them to orderlist
                auxiliary_dict["orderlist_series"] = pd.concat([auxiliary_dict["orderlist_series"], unordered_whitelist_series], ignore_index=True)

        for df in auxiliary_dict["len_df"], auxiliary_dict["centromere_df"], auxiliary_dict["highlight_df"]:
            if df is not None:
                df.rename(index=auxiliary_dict["syn_dict"], inplace=True)
        auxiliary_dict["orderlist_series"] = auxiliary_dict["orderlist_series"][::-1]
        return tmp_df


class NewlinePreservingArgParserHelpFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        # Split on explicit newlines, but wrap each paragraph normally
        lines = text.splitlines()
        wrapped = [textwrap.fill(line, width, initial_indent=indent, subsequent_indent=indent) for line in lines]
        return "\n".join(wrapped)

    def _split_lines(self, text, width):
        lines = text.splitlines()
        wrapped_lines = []
        for line in lines:
            wrapped_lines.extend(textwrap.wrap(line, width))
            if not line:
                wrapped_lines.append("")
        return wrapped_lines + [""]
