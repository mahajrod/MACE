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

    @staticmethod
    def get_filenames_for_extension(dir_path, extension_list, force_uniq=True):
        filelist = []
        for extension in extension_list:
            filelist += list(glob.glob(str(dir_path) + "/*{0}".format(extension)))
        if not filelist:
            return None
        # print(filelist)
        if force_uniq:
            if len(filelist) > 1:
                raise ValueError("Found more than one file with extensions: {0} in directory {1}".format(",".join(extension_list),
                                                                                                         str(dir_path)))
            else:
                return filelist[0]

        return filelist

    # ---------- Methods for parsing config file for draw_macrosynteny.py -------------
    @staticmethod
    def check_for_inversion_and_strand_switch_request(scaffold_id, strand_switch_symbol, inversion_symbol, ):
        request_dict = {"inversion": False,
                        "query_switch_strand": False,
                        "target_switch_strand": False}
        processed_scaffold_id = scaffold_id
        if scaffold_id[-1] == inversion_symbol:
            request_dict["inversion"] = True
            processed_scaffold_id = scaffold_id[:-1]
        if processed_scaffold_id[-3:] == strand_switch_symbol * 3:
            request_dict["query_switch_strand"] = True
            request_dict["target_switch_strand"] = True
            processed_scaffold_id = processed_scaffold_id[:-3]
        elif processed_scaffold_id[-2:] == strand_switch_symbol * 2:
            request_dict["target_switch_strand"] = True
            processed_scaffold_id = processed_scaffold_id[:-2]
        elif processed_scaffold_id[-1:] == strand_switch_symbol * 1:
            request_dict["query_switch_strand"] = True
            processed_scaffold_id = processed_scaffold_id[:-1]

        if len(processed_scaffold_id) == 0:
            raise ValueError(f"ERROR!!! No symbols were left after processing of the scaffold id {scaffold_id}!")

        return processed_scaffold_id, request_dict

    @staticmethod
    def get_original_scaffold_id(syn_df, scaffold_id):
        if scaffold_id not in list(syn_df["syn"]):
            # print(syn_df["syn"])
            return scaffold_id
        return syn_df.index[list(syn_df["syn"]).index(scaffold_id)]
    # ---------- End of Methods for parsing config file for draw_macrosynteny.py -------------

    def read_mace_auxiliary_input(self,
                                  len_file=None,
                                  whitelist_file=None, max_scaffolds=50,
                                  orderlist_file=None,
                                  invertlist_file=None, inverted_scaffold_label="'",
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
                                  median_coverage=None, mean_coverage=None,
                                  scaffold_color_file=None,
                                  ):
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
        auxiliary_dict["preinvert_len_df"] = deepcopy(auxiliary_dict["len_df"])

        if legend_file is not None:
            if os.path.exists(legend_file):
                try:  # legend_file might be empty
                    auxiliary_dict["legend_df"] = pd.read_csv(legend_file, header=None, index_col=0, usecols=[0, 1],
                                                              names=["entry", "color"], dtype={"entry": str, "color": str},
                                                              sep="\t", comment=None)
                    #len_df.index = pd.Index(list(map(str, len_df.index)))
                except pd.errors.EmptyDataError:
                    auxiliary_dict["legend_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Legend file {legend_file} doesn't exist!")
        else:
            auxiliary_dict["legend_df"] = None

        if scaffold_color_file is not None:
            if os.path.exists(scaffold_color_file):
                try:  # legend_file might be empty
                    auxiliary_dict["scaffold_color_df"] = pd.read_csv(scaffold_color_file, header=None, index_col=0, usecols=[0, 1],
                                                                      names=["scaffold", "color"],
                                                                      dtype={"scaffold": str, "color": str},
                                                                      sep="\t", comment=None)
                    #len_df.index = pd.Index(list(map(str, len_df.index)))
                except pd.errors.EmptyDataError:
                    auxiliary_dict["scaffold_color_df"] = None
            else:
                raise FileNotFoundError(f"ERROR!!! Scaffold color file {scaffold_color_file} doesn't exist!")
        else:
            auxiliary_dict["scaffold_color_df"] = None

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
        auxiliary_dict["inverted_scaffold_label"] = inverted_scaffold_label

        return auxiliary_dict

    def update_genome_dict_from_genome_config(self, genome_auxiliary_dict, genome_config,
                                              strand_switch_label="*", inverted_scaffold_label="'",):
        # genome_auxiliary_dict has following definition
        # genome_auxiliary_dict = {"genomes": OrderedDict(),
        #                          "genome_colors": OrderedDict()}
        #genome_auxiliary_dict["genome_color_df"] = {"genome_id": [],
        #                                            "color": []}
        if genome_config is not None:
            genome_config_tmp_dict = OrderedDict()
            # read genome colors first
            genome_auxiliary_dict["genome_color_df"] = pd.read_csv(genome_config, sep="\t", usecols=[0, 1], na_values=".",
                                                                   names=["genome", "color"], index_col=["genome"])

            # read scaffold_settings
            with open(genome_config, "r") as in_fd:
                for line in in_fd:
                    if line == "\n":
                        continue
                    line_list = line.strip().split("\t", 2)
                    #genome_auxiliary_dict["genome_color_df"]["genome_id"].append(line_list[0])
                    #genome_auxiliary_dict["genome_color_df"]["color"].append(line_list[1] if line_list[1] != "." else None)  # read_colors
                    genome_config_tmp_dict[line_list[0]] = list(map(lambda s: s.split(","), line_list[2].split("\t")))
            #genome_auxiliary_dict["genome_color_df"] = pd.DataFrame.from_dict(genome_auxiliary_dict["genome_color_df"]).set_index("genome_id")

            for genome in genome_auxiliary_dict["genomes"]:
                for series_entry in "invertlist_series", "orderlist_series", "queryswitchstrandlist_series", "targetswitchstrandlist_series":
                    genome_auxiliary_dict["genomes"][genome][series_entry] = []

                for entry in genome_config_tmp_dict[genome]:
                    for scaffold_id in entry:
                        processed_scaffold_id, request_dict = self.check_for_inversion_and_strand_switch_request(scaffold_id,
                                                                                                                 strand_switch_label,
                                                                                                                 inverted_scaffold_label)
                        genome_auxiliary_dict["genomes"][genome]["orderlist_series"].append(processed_scaffold_id)
                        if request_dict["inversion"]:
                            genome_auxiliary_dict["genomes"][genome]["invertlist_series"].append(processed_scaffold_id)

                        if request_dict["query_switch_strand"]:
                            genome_auxiliary_dict["genomes"][genome]["queryswitchstrandlist_series"].append(self.get_original_scaffold_id(genome_auxiliary_dict["genomes"][genome]["syn_dict"],
                                                                                                            processed_scaffold_id))
                        if request_dict["target_switch_strand"]:
                            genome_auxiliary_dict["genomes"][genome]["targetswitchstrandlist_series"].append(self.get_original_scaffold_id(genome_auxiliary_dict["genomes"][genome]["syn_dict"],
                                                                                                             processed_scaffold_id))

                genome_auxiliary_dict["genomes"][genome]["invertlist_series"] = pd.Series(genome_auxiliary_dict["genomes"][genome]["invertlist_series"], dtype='str')
                genome_auxiliary_dict["genomes"][genome]["orderlist_series"] = pd.Series(genome_auxiliary_dict["genomes"][genome]["orderlist_series"], dtype='str')
                genome_auxiliary_dict["genomes"][genome]["queryswitchstrandlist_series"] = pd.Series(genome_auxiliary_dict["genomes"][genome]["queryswitchstrandlist_series"],
                                                                                                     dtype='str')
                genome_auxiliary_dict["genomes"][genome]["targetswitchstrandlist_series"] = pd.Series(genome_auxiliary_dict["genomes"][genome]["targetswitchstrandlist_series"],
                                                                                                      dtype='str')

    def resolve_mace_single_genome_input(self, auxiliary_dict, records_df=None):
        if auxiliary_dict["len_df"] is None:
            raise ValueError("ERROR!!! Length dataframe was not parsed!!!")

        if records_df is None:
            tmp_df = records_df
        else:
            if (auxiliary_dict["whitelist_series"] is not None) and auxiliary_dict["whitelist_series"].empty:
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

        if (auxiliary_dict["whitelist_series"] is not None) and (not auxiliary_dict["whitelist_series"].empty):
            # remove from orderlist scaffolds which are absent in the whitelist
            if (auxiliary_dict["orderlist_series"] is not None) and auxiliary_dict["orderlist_series"].empty:
                auxiliary_dict["orderlist_series"] = auxiliary_dict["whitelist_series"].replace(auxiliary_dict["syn_dict"])
            else:
                renamed_whitelist_series = auxiliary_dict["whitelist_series"].replace(auxiliary_dict["syn_dict"])
                # remove orderlist scaffolds that are not in whitelist
                auxiliary_dict["orderlist_series"] = auxiliary_dict["orderlist_series"].loc[auxiliary_dict["orderlist_series"].isin(renamed_whitelist_series)]
                # find whitelist scaffolds that are not in orderlist
                unordered_whitelist_series = renamed_whitelist_series[~renamed_whitelist_series.isin(auxiliary_dict["orderlist_series"])]
                # add them to orderlist
                auxiliary_dict["orderlist_series"] = pd.concat([auxiliary_dict["orderlist_series"], unordered_whitelist_series], axis="index", ignore_index=True)
        if auxiliary_dict["orderlist_series"] is not None:
            auxiliary_dict["orderlist_series"] = auxiliary_dict["orderlist_series"][::-1]

        if (auxiliary_dict["scaffold_color_df"] is not None) and (not auxiliary_dict["scaffold_color_df"].empty):
            reordered_color_df_list = []
            color_in_orderlist_df = auxiliary_dict["scaffold_color_df"].loc[auxiliary_dict["orderlist_series"]]
            color_not_in_orderlist_df = auxiliary_dict["scaffold_color_df"].loc[~auxiliary_dict["scaffold_color_df"].index.isin(auxiliary_dict["orderlist_series"])]
            if not color_in_orderlist_df.empty:
                reordered_color_df_list.append(color_in_orderlist_df)
            if not color_not_in_orderlist_df.empty:
                reordered_color_df_list.append(color_not_in_orderlist_df)
            auxiliary_dict["scaffold_color_df"] = pd.concat(reordered_color_df_list, axis="index")

        for df in auxiliary_dict["len_df"], auxiliary_dict["centromere_df"], auxiliary_dict["highlight_df"]:
            if df is not None:
                df.rename(index=auxiliary_dict["syn_dict"], inplace=True)

        auxiliary_dict["preinvert_len_df"] = deepcopy(auxiliary_dict["len_df"])
        
        if auxiliary_dict["invertlist_series"] is not None:
            if auxiliary_dict["len_df"] is not None:
                auxiliary_dict["len_df"] = auxiliary_dict["len_df"].rename(index=dict(zip(auxiliary_dict["invertlist_series"],
                                                                                          [scaf + auxiliary_dict["inverted_scaffold_label"] for scaf in auxiliary_dict["invertlist_series"]])))

            if auxiliary_dict["orderlist_series"] is not None:
                auxiliary_dict["orderlist_series"] = auxiliary_dict["orderlist_series"].replace(dict(zip(auxiliary_dict["invertlist_series"],
                                                                                                     [scaf + auxiliary_dict["inverted_scaffold_label"] for scaf in auxiliary_dict["invertlist_series"]])))
            if auxiliary_dict["scaffold_color_df"] is not None:
                auxiliary_dict["scaffold_color_df"] = auxiliary_dict["scaffold_color_df"].rename(index=dict(zip(auxiliary_dict["invertlist_series"],
                                                                                                                [scaf + auxiliary_dict["inverted_scaffold_label"] for scaf in auxiliary_dict["invertlist_series"]])))
            if auxiliary_dict["centromere_df"] is not None:
                auxiliary_dict["centromere_df"] = self.invert_coordinates_in_region_table(auxiliary_dict["centromere_df"],
                                                                                          auxiliary_dict["invertlist_series"],
                                                                                          auxiliary_dict["len_df"],
                                                                                          "scaffold", "start", "end",
                                                                                          inverted_scaffolds_label=auxiliary_dict["inverted_scaffold_label"])
        
        return tmp_df

    @staticmethod
    def invert_coordinates_in_synteny_table(df, scaffold_list, length_df, scaffold_column, start_column, end_column,
                                            strand_column, inverted_scaffolds_label="'"):
        temp_df = deepcopy(df)
        columns_list = list(temp_df.columns)
        temp_df["length_column"] = 0
        original_index_name = temp_df.index.name

        temp_df.reset_index(drop=False, inplace=True)
        temp_df.set_index(scaffold_column, inplace=True)
        # temp_df.to_csv("tmp", sep="\t", index=True, header=True)
        for scaffold in temp_df.index.unique():
            temp_df.loc[scaffold, "length_column"] = length_df.loc[scaffold, "length"]

        (temp_df.loc[temp_df.index.isin(scaffold_list), start_column],
         temp_df.loc[temp_df.index.isin(scaffold_list), end_column]) = temp_df.loc[temp_df.index.isin(scaffold_list), "length_column"] - temp_df.loc[temp_df.index.isin(scaffold_list), end_column], \
                                                                       temp_df.loc[temp_df.index.isin(scaffold_list), "length_column"] - temp_df.loc[temp_df.index.isin(scaffold_list), start_column]

        plus_indexes, minus_indexes = (temp_df[strand_column] == "+") & temp_df.index.isin(scaffold_list), (
                    temp_df[strand_column] == "-") & temp_df.index.isin(scaffold_list)
        temp_df.loc[plus_indexes, strand_column], temp_df.loc[minus_indexes, strand_column] = "-", "+"
        temp_df.reset_index(drop=False, inplace=True)
        if inverted_scaffolds_label is not None:
            for scaffold in scaffold_list:
                temp_df.loc[temp_df[scaffold_column] == scaffold, scaffold_column] = scaffold + inverted_scaffolds_label
        if (original_index_name is not None) and original_index_name in temp_df.columns:
            temp_df.set_index(original_index_name, inplace=True)

        return temp_df[columns_list]  # remove added length column and restore column order

    """# function from draw_synteny.py, slightly differs from above method 
    def invert_coordinates_in_synteny_table(df, scaffold_list, length_df, scaffold_column, start_column, end_column, strand_column, inverted_scaffolds_label="'"):
        temp_df = deepcopy(df)
        columns_list = list(temp_df.columns)
        temp_df["length_column"] = 0
        temp_df.set_index(scaffold_column, inplace=True)

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
        return temp_df[columns_list]  # remove added length column and restore column order
    """

    @staticmethod
    def invert_coordinates_in_region_table(df, scaffold_list, length_df, scaffold_column, start_column, end_column,
                                           inverted_scaffolds_label="'"):
        if df.empty:
            return df
        temp_df = deepcopy(df)

        if temp_df.index.name != scaffold_column:
            temp_df.reset_index(inplace=True, drop=False)
            temp_df.set_index(scaffold_column, inplace=True)

        columns_list = list(temp_df.columns)
        for scaffold in temp_df.index.unique():
            if scaffold in length_df.index:
                temp_df.loc[scaffold, "length_column"] = length_df.loc[scaffold, "length"]

        temp_df.loc[temp_df.index.isin(scaffold_list), start_column], temp_df.loc[
            temp_df.index.isin(scaffold_list), end_column] = temp_df.loc[temp_df.index.isin(scaffold_list), "length_column"] - temp_df.loc[temp_df.index.isin(scaffold_list), end_column], \
                                                             temp_df.loc[temp_df.index.isin(scaffold_list), "length_column"] - temp_df.loc[temp_df.index.isin(scaffold_list), start_column]

        temp_df.reset_index(drop=False, inplace=True)
        if inverted_scaffolds_label is not None:
            for scaffold in scaffold_list:
                temp_df.loc[temp_df[scaffold_column] == scaffold, scaffold_column] = scaffold + inverted_scaffolds_label

        temp_df.set_index(scaffold_column, inplace=True)
        return temp_df[columns_list]  # remove added length column and restore column order

    @staticmethod
    def expand_path(path_template: str, skip=False):
        if (path_template[0] == "/") or skip:
            print("Skipping expanding path {0} as it is global path or expanding was turned off ...".format(
                path_template))
            # avoid expanding global paths
            return path_template

        path_list = list(Path("./").glob(path_template))
        if len(path_list) > 1:
            raise ValueError(
                "ERROR!!! There is more than one file corresponding to the template {0} ...".format(path_template))
        elif len(path_list) == 0:
            raise ValueError("ERROR!!! There are no files corresponding to the template {0} ...".format(path_template))
        else:
            return str(path_list[0])


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

