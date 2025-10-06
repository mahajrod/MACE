#!/usr/bin/env python
__author__ = 'mahajrod'

import re
import os
import sys
from copy import deepcopy

from collections.abc import Iterable
from collections import OrderedDict

import pandas as pd
import matplotlib.colors as mpl_colors

if sys.version_info[0] == 3:
    from io import TextIOWrapper as file

implemented_tracktype_dict = {
                              "marker":   [],   # track to this dictionary to register it
                              "plot": [],
                              "window":   [],
                              "hist": [],
                              }

recognizable_named_color_list = set(mpl_colors.TABLEAU_COLORS) | set(mpl_colors.CSS4_COLORS) | set(mpl_colors.XKCD_COLORS) | set(["default", ])


full_marker_type_list = ["rectangle", "ellipse", "circle"]
short_marker_type_list = ["r", "e", "c"]
recognizable_marker_type_list = full_marker_type_list + short_marker_type_list


class DataDF(pd.DataFrame):
    _metadata = ["START_COL_INDEX", "END_COL_INDEX", "parameter_separator", "subtype_separator",
                 "attachment_separator", "default_tracktype", "check_method_dict", "convert_method_dict",
                 "track_number", "track_list", "autodetected_tracktype_list", "defaulted_tracktype_list",
                 "track_type_dict", "track_parameter_dict", "track_attachment_dict", "track_metadata_df", "parsed",
                 "attached_track_dict", "attachment_relation_dict", "track_feature_dict",
                 "track_name_list", "track_type_list"]

    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=True,
                 syn_df=None, whitelist_series=None, blacklist_series=None,
                 empty=False, parsed=False, default_tracktype=None, track_name_list=None, track_type_list=None,
                 parameter_separator="&", subtype_separator="$", attachment_separator="@"):
        
        pd.DataFrame.__init__(self, data=data, index=index, columns=columns, dtype=dtype, copy=copy,)
        #self = pd.DataFrame(data=data, index=index, columns=columns, dtype=dtype, copy=copy)

        #---- constants ----
        self.START_COL_INDEX = 0
        self.END_COL_INDEX = 1

        #---- parameters and settings ----
        self.parameter_separator = parameter_separator
        self.subtype_separator = subtype_separator
        self.attachment_separator = attachment_separator
        self.default_tracktype = default_tracktype

        self.check_method_dict = {  # for a new track type code a check method and add it here
                                  "marker":     self.check_if_marker,
                                  "plot":       self.check_if_plot,
                                  "window":     self.check_if_window,
                                  "hist":       self.check_if_hist,
                                  }

        self.convert_method_dict = {  # for a new attribute or track type code a check method and add it here, parsing function must return a dataframe with two level Multiindex as column names
                                  "marker":     self.parse_string_column,
                                  "plot":       self.parse_float_column,
                                  "window":     self.parse_string_column,
                                  "hist":       self.parse_float_list_column,
                                  "colors":     self.parse_string_list_column,
                                  "bg_color":   self.parse_string_column,
                                  }

        #---- default metadata (is filled during parsing) ----
        self.track_number = None
        self.track_list = []
        self.track_name_list = track_name_list if track_name_list is not None else []
        self.track_type_list = track_type_list if track_type_list is not None else []
        self.autodetected_tracktype_list = []
        self.defaulted_tracktype_list = []
        self.track_type_dict = {}
        self.track_parameter_dict = {}
        self.track_attachment_dict = {}
        self.track_metadata_df = None
        self.parsed = parsed
        self.attached_track_dict = {}
        self.attachment_relation_dict = {}
        self.track_feature_dict = {}

        if (whitelist_series is not None) or (blacklist_series is not None):
            #print("AAAAAA")
            self.filter_scaffolds(whitelist_series=whitelist_series,
                                  blacklist_series=blacklist_series,
                                  inplace=True)

        if (syn_df is not None) and (not syn_df.empty):
            self.rename_scaffolds(self, syn_df)

    @property
    def _constructor(self):
        return DataDF

    def parse(self):
        self._check_input_data(empty=self.empty)
        self.parse_metadata_from_column_names()
        self.parse_features()
        self.parse_data()
        self.parsed = True

        return self

    def _check_input_data(self, empty=False):
        #print(column_names)
        # input data should have scaffold_ids as index, and two first columns should be start and stop of the region/feature, i.e. integers

        # check if index of data is of string type
        if not pd.api.types.is_string_dtype(self.index):
            raise ValueError("ERROR!!! Index of DataDF, i.e. scaffold ids, of the input data contains non string values!")

        # check that data have three columns (two in case of empty track)
        #print(self.columns)
        if empty:
            if len(self.columns) < 2:
                raise ValueError("ERROR!!! Empty DataDF must contain at least two columns (start and stop)")
        else:
            if len(self.columns) < 3:
                print(self)
                raise ValueError("ERROR!!! Non empty DataDF must contain at least three columns (start, stop and values for at least a single track)")
        # check if start (first) and stop (second) columns contain only integer values
        if not pd.api.types.is_integer_dtype(self.iloc[:, 0]):
            raise ValueError("ERROR!!! First column of DataDF, i.e. 'start', of the input data contains non integer values!")

        if not pd.api.types.is_integer_dtype(self.iloc[:, 1]):
            raise ValueError("ERROR!!! Second column of DataDF, i.e. 'end', of the input data contains non integer values!")

        # check column names of track values and parameters
        if not empty:
            for column_name in self.columns[2:]:
                if column_name.count(self.parameter_separator) > 1:
                    raise ValueError(f"ERROR!!! At least one column name contain more than one parameter separator '{self.parameter_separator}'. "
                                     f"Change parameter separator or rename corresponding column.")
                if column_name.count(self.attachment_separator) > 1:
                    raise ValueError(f"ERROR!!! At least one column name contain more than one attachment separator '{self.attachment_separator}'. "
                                     f"Change attachment separator or rename corresponding column.")

        if self.default_tracktype is not None:
            if not (self.default_tracktype in implemented_tracktype_dict):
                raise ValueError(f"ERROR!!! Unrecognized default track type ({self.default_tracktype}). "
                                 f"Allowed track types: {','.join(list(implemented_tracktype_dict.keys()))}. ")

    def parse_metadata_from_column_names(self):  # method assumes that DataDF is not empty

        column_names = self.columns[2:]

        header_parsing_list = []  # [ <track_name>, <track_parameter>, <parameter_value>, <feature_list>, <attachment>]

        # Parsing of the metadata from column names starts from the end of the column name
        #print(column_names)
        for column_name in column_names.get_level_values(0) if isinstance(column_names, pd.MultiIndex) else column_names:
            #print(column_name)
            # extracting_attachement
            attachement_parsing_list = column_name.split(self.attachment_separator)
            if len(attachement_parsing_list) == 1:
                attachment = pd.NA
            else:
                attachment = attachement_parsing_list[-1]
            column_name_prefix = attachement_parsing_list[0]

            # extracting subtypes/features
            feature_parsing_list = column_name_prefix.split(self.subtype_separator)
            if len(feature_parsing_list) == 1:
                feature_list = pd.NA
            else:
                feature_list = feature_parsing_list[1:]

            column_name_pre_prefix = feature_parsing_list[0]

            # extracting type/parameter
            parameter_parsing_list = column_name_pre_prefix.split(self.parameter_separator)
            if len(parameter_parsing_list) == 1:
                track_parameter = "type"
                parameter_value = pd.NA
            else:
                if parameter_parsing_list[1] in implemented_tracktype_dict:
                    track_parameter = "type"
                    parameter_value = parameter_parsing_list[1]
                elif parameter_parsing_list[1] == "default":
                    track_parameter = "type"
                    parameter_value = "default"
                else:
                    track_parameter = "parameter"
                    parameter_value = parameter_parsing_list[1]
            track_name = parameter_parsing_list[0]

            header_parsing_list.append([track_name, track_parameter, parameter_value, feature_list, attachment, column_name])

        parsing_df = pd.DataFrame.from_records(header_parsing_list,
                                               columns=["track_name", "parameter", "parameter_value", "feature_list",
                                                        "attachment", "column_name"],
                                               )

        parsing_df.set_index("track_name", inplace=True)
        type_df = parsing_df[parsing_df["parameter"] == "type"]
        absent_type_df = type_df[type_df["parameter_value"].isna()]

        if not absent_type_df.empty:
            sys.stderr.write(f"WARNING!!! Some tracks({', '.join(absent_type_df.index)}) don't have type preset. "
                             f"Trying auto detection.\n")
            for trackname in absent_type_df.index:
                sys.stderr.write(f"{trackname}:\n")
                detected_tracktype = self.detect_track_type(absent_type_df.loc[trackname, "column_name"])
                if detected_tracktype is not None:
                    self.autodetected_tracktype_list.append(trackname)
                    type_df.loc[trackname, "parameter_value"] = detected_tracktype
                else:
                    type_df.loc[trackname, "parameter_value"] = "default"

        default_type_df = type_df[type_df["parameter_value"] == "default"]
        if not default_type_df.empty:
            sys.stderr.write(f"WARNING!!! Some tracks ({', '.join(default_type_df.index)}) have preset type 'default' or failed auto detection."
                             f"Assigning default track type to them...\n")
            #print(default_type_df)
            if self.default_tracktype is None:
                raise ValueError(f"ERROR!!! Default track type is not set. Set types for corresponding tracks ({', '.join(default_type_df.index)}) manually in the header")
            type_df.loc[default_type_df.index, "parameter_value"] = self.default_tracktype
            self.defaulted_tracktype_list += list(default_type_df.index)

        parsing_df.loc[parsing_df["parameter"] == "type", "parameter_value"] = type_df["parameter_value"]
        #print(type_df)
        sys.stdout.write("Results of track metadata parsing:\n")
        sys.stdout.write(str(parsing_df) + "\n")

        self.track_name_list = list(parsing_df[parsing_df["parameter"] == "type"].index)
        self.track_type_list = list(parsing_df[parsing_df["parameter"] == "type"]["parameter_value"])
        self.track_metadata_df = parsing_df
        return parsing_df

    def detect_track_type(self, track_column_name):
        for track_type in self.check_method_dict:
            sys.stderr.write(f"\tIs {track_type}?:\t")
            if self.check_method_dict[track_type](track_column_name):
                sys.stderr.write("True\n")
                return track_type
            sys.stderr.write("False\n")
        return None
        #else:
        #    raise ValueError(f"ERROR!!! Track type was unrecognized for {track_column_name}! Either set it directly in the column name or add support for this data type to the code.")

    def check_if_hist(self, track_column_name):
        # values in column should be strings containing comma-separated list of float values
        is_hist = False
        if not pd.api.types.is_string_dtype(self[track_column_name].dropna()):
            return False
        # try to convert first nonNA value in the column
        value_list = self[track_column_name].dropna().iloc[0].split(",")
        if len(value_list) == 1:
            return False
        for converter in [float]:
            try:
                map(converter, value_list)
            except:
                is_hist = is_hist | False
            else:
                is_hist = is_hist | True

        return is_hist

    def check_if_plot(self, track_column_name):
        # values in column should be numerical or strings convertable to a single numerical value
        is_plot = False
        if not (pd.api.types.is_numeric_dtype(self[track_column_name].dropna()) or pd.api.types.is_string_dtype(self[track_column_name].dropna())):
            return False

        # try to convert first non NA value in column

        for converter in float,:
            try:
                converter(self[track_column_name].dropna().iloc[0])
            except:
                is_plot = is_plot | False
            else:
                is_plot = is_plot | True

        return is_plot

    def check_if_window(self, track_column_name):
        # values in column should encode a single color. It maybe a color name recognizable by matplotlib or a six or eight letter (RGB or RGBA) hex color code preceding with #
        is_window = False
        if not pd.api.types.is_string_dtype(self[track_column_name].dropna()):
            return False
        # try to convert first non NA value
        if self[track_column_name].dropna().iloc[0] in recognizable_named_color_list:
            return True

        if "," in self[track_column_name].dropna().iloc[0]:  # to avoid mess with hist
            return False

        for converter in mpl_colors.to_rgb, mpl_colors.to_rgba:
            try:
                list(map(converter, self[track_column_name].dropna().iloc[0].split(",")))
            except:
                is_window = is_window | False
            else:
                #print(map(converter, self[track_column_name].iloc[0].split(",")))
                is_window = is_window | True

        return is_window

    def check_if_marker(self, track_column_name):
        # values in column should encode a recognizable marker type:
        if self[track_column_name].dropna().iloc[0] in recognizable_marker_type_list:
            #print(self[track_column_name].iloc[0])
            return True
        return False

    def parse_float_column(self, track_column_name):
        tmp_df = self[[track_column_name]].copy().astype(float)
        tmp_df.columns = pd.MultiIndex.from_tuples([(track_column_name, track_column_name)],
                                                   names=["parameter", "element"])
        return tmp_df

    def parse_int_column(self, track_column_name):    # using nullable Int64 type
        tmp_df = self[[track_column_name]].copy().astype('Int64')
        tmp_df.columns = pd.MultiIndex.from_tuples([(track_column_name, track_column_name)],
                                                   names=["parameter", "element"])
        return tmp_df

    def parse_string_column(self, track_column_name):
        tmp_df = self[[track_column_name]].copy()
        tmp_df.columns = pd.MultiIndex.from_tuples([(track_column_name, track_column_name)],
                                                   names=["parameter", "element"])
        return tmp_df

    def parse_float_list_column(self, track_column_name):
        tmp_df = self[track_column_name].apply(lambda s: pd.Series(map(lambda f: pd.NA if self.check_na(f) else float(f),
                                                                       [pd.NA] if pd.isna(s) else s.split(","))))
        tmp_df.columns = pd.MultiIndex.from_product([[track_column_name], list(map(lambda s: f'c{s}', list(tmp_df.columns)))], names=["parameter", "element"])
        return tmp_df

    def parse_int_list_column(self, track_column_name):             # using nullable Int64 type
        tmp_df = self[track_column_name].apply(lambda s: pd.Series(map(lambda f: pd.NA if self.check_na(f) else float(f),
                                                                   [pd.NA] if pd.isna(s) else s.split(",")))).astype('Int64')
        tmp_df.columns = pd.MultiIndex.from_product([[track_column_name], list(map(lambda s: f'c{s}', list(tmp_df.columns)))], names=["parameter", "element"])
        return tmp_df

    def parse_string_list_column(self, track_column_name):
        tmp_df = self[track_column_name].apply(lambda s: pd.Series(pd.NA if pd.isna(s) else s.split(",")))
        tmp_df.columns = pd.MultiIndex.from_product([[track_column_name], list(map(lambda s: f'c{s}', list(tmp_df.columns)))], names=["parameter", "element"])

        return tmp_df

    @staticmethod
    def check_na(value):
        na_check = pd.isnull(value)
        return sum(na_check) if isinstance(na_check, Iterable) else na_check

    @staticmethod
    def get_filtered_entry_list(entry_list,
                                entry_black_list=None,
                                sort_entries=False,
                                entry_ordered_list=None,
                                entry_white_list=None,
                                invert_match=False):
        white_set = set(entry_white_list) if entry_white_list is not None else set()
        black_set = set(entry_black_list) if entry_black_list is not None else set()
        entry_set = set(entry_list)

        if white_set:
            filtered_entry_set = entry_set & white_set
        if black_set:
            filtered_entry_set = entry_set - black_set
        if (not white_set) and (not black_set):
            filtered_entry_set = set()

        if invert_match:
            filtered_entry_set = entry_set - filtered_entry_set
            return list(filtered_entry_set)
        else:
            filtered_entry_list = list(filtered_entry_set)
            if sort_entries:
                filtered_entry_list.sort()

            final_entry_list = []

            if entry_ordered_list:
                for entry in entry_ordered_list:
                    if entry in filtered_entry_list:
                        final_entry_list.append(entry)
                        filtered_entry_list.remove(entry)
                    else:
                        print("WARNING!!!Entry(%s) from order list is absent in list of entries!" % entry)
                return final_entry_list + filtered_entry_list
            else:
                return filtered_entry_list

    def filter_scaffolds(self, whitelist_series=None, blacklist_series=None, inplace=False):
        if inplace:
            filtered_scaffolds = self.get_filtered_entry_list(self.index,
                                                              entry_black_list=blacklist_series,
                                                              entry_white_list=whitelist_series,
                                                              invert_match=True)
            self.drop(filtered_scaffolds, inplace=True)
        else:
            filtered_scaffolds = self.get_filtered_entry_list(self.index,
                                                              entry_black_list=blacklist_series,
                                                              entry_white_list=whitelist_series,
                                                              invert_match=False)
            return self.loc[self.index.isin(filtered_scaffolds)]

    def rename_scaffolds(self, syn_df, syn_column="syn"):
        self.rename(index=syn_df[syn_column].to_dict(), inplace=True)

    def parse_data(self):
        column_names = list(self.columns)
        columns_to_drop = column_names[2:]
        column_multiindex = pd.MultiIndex.from_arrays([["start", "end"], ["start", "end"]], names=["parameter", "element"])

        new_column_index = 0
        for row_index in range(0, len(self.track_metadata_df)):
            parameter_value = self.track_metadata_df["parameter_value"].iloc[row_index]
            column_name = self.track_metadata_df["column_name"].iloc[row_index]

            tmp = self.convert_method_dict[parameter_value](column_name)

            column_multiindex = column_multiindex.append(tmp.columns)
            tmp.columns = tmp.columns.droplevel(level=0)
            #print(tmp)
            for element in tmp.columns:
                self[f"new{new_column_index}"] = tmp[element]
                new_column_index += 1
        #print(self)
        self.drop(columns_to_drop, axis=1, inplace=True)
        #print(column_multiindex)
        self.columns = column_multiindex

    def parse_features(self):
        for track_name in list(self.track_metadata_df[(self.track_metadata_df["parameter"] == "type")].index):
            self.track_feature_dict[track_name] = OrderedDict()
            for parameter, parameter_value, feature_list, column_name in self.track_metadata_df.loc[[track_name],
                                                                                                    ["parameter",
                                                                                                     "parameter_value",
                                                                                                     "feature_list",
                                                                                                     "column_name"]].itertuples(index=False):
                tmp_list = list(map(lambda s: s.split("="), feature_list)) if not self.check_na(feature_list) else None
                feature_group = "general" if parameter == "type" else parameter_value
                self.track_feature_dict[track_name][feature_group] = OrderedDict()
                if tmp_list is not None:
                    for entry in tmp_list:
                        if len(entry) == 1:
                            self.track_feature_dict[track_name][feature_group][entry[0]] = True
                        elif len(entry) == 2:
                            try:  # try converting feature value first to int, and after that to float
                                entry[1] = int(entry[1])
                            except:
                                try:
                                    entry[1] = float(entry[1])
                                except:
                                    pass
                            self.track_feature_dict[track_name][feature_group][entry[0]] = entry[1]

                        else:
                            raise ValueError(f"ERROR!!! More than one feature separators('=') in one of the columns")



    """    
    def get_host_chain_for_track(self, track_name, host_chain_list=[]):
        host_track_name = self.track_metadata_df[self.track_metadata_df["parameter"] == "type"].loc[track_name, "attachment"]
        if self.check_na(host_track_name):
            return (host_chain_list + [track_name])[::-1]
        else:
            return self.get_host_chain_for_track(host_track_name, host_chain_list + [track_name])

    def update_nested_dict_from_list(self, chain_list, nested_dict):
        print(chain_list)
        if chain_list[0] not in nested_dict:
            nested_dict[chain_list[0]] = OrderedDict()
        if len(chain_list) == 1:
            return 0
        return self.update_nested_dict_from_list(chain_list[1:], nested_dict[chain_list[0]])

    def resolve_attachment_relations(self):
        tmp_dict = OrderedDict()
        for track_name in list(self.track_metadata_df.index.unique()):
            host_chain_list = self.get_host_chain_for_track(track_name)
            print(track_name, host_chain_list)
            self.update_nested_dict_from_list(host_chain_list, tmp_dict)
        self.attachment_relation_dict = tmp_dict
        return tmp_dict
    """
    def split_tracks(self, resolve_attachments=True):
        track_list = []

        if resolve_attachments:
            independent_tracks = list(self.track_metadata_df[(self.track_metadata_df["parameter"] == "type") & (self.track_metadata_df["attachment"].isna())].index)
            for track_name in independent_tracks:
                track_list.append(DataDF(self[["start", "end"] + list(self.track_metadata_df.loc[[track_name]]["column_name"])], parsed=True))
                track_list[-1].track_metadata_df = self.track_metadata_df.loc[[track_name]].copy()
                track_list[-1].track_feature_dict = {track_name: deepcopy(self.track_feature_dict[track_name])}
                track_list[-1].track_name_list = [track_name]
                track_list[-1].track_type_list = [self.track_metadata_df[(self.track_metadata_df["parameter"] == "type")].loc[track_name, "parameter_value"]]
                for attached_track_name in list(self.track_metadata_df[(self.track_metadata_df["parameter"] == "type") & (self.track_metadata_df["attachment"] == track_name)].index):
                    track_list[-1].attached_track_dict[attached_track_name] = DataDF(self[["start", "end"] + list(self.track_metadata_df.loc[[attached_track_name]]["column_name"])], parsed=True)
                    track_list[-1].attached_track_dict[attached_track_name].track_metadata_df = self.track_metadata_df.loc[[attached_track_name]].copy()
                    track_list[-1].attached_track_dict[attached_track_name].track_feature_dict = {attached_track_name: deepcopy(self.track_feature_dict[attached_track_name])}
                    track_list[-1].attached_track_dict[attached_track_name].track_name_list = [attached_track_name]
                    track_list[-1].attached_track_dict[attached_track_name].track_type_list = [self.track_metadata_df[(self.track_metadata_df["parameter"] == "type")].loc[attached_track_name, "parameter_value"]]
        else:
            for track_name in list(self.track_metadata_df.index.unique()):
                track_list.append(DataDF(self[["start", "end"] + list(self.track_metadata_df.loc[[track_name]]["column_name"])], parsed=True))
                track_list[-1].track_metadata_df = self.track_metadata_df.loc[[track_name]].copy()
                track_list[-1].track_feature_dict = {track_name: deepcopy(self.track_feature_dict[track_name])}
                track_list[-1].track_name_list = [track_name]
                track_list[-1].track_type_list = [self.track_metadata_df[(self.track_metadata_df["parameter"] == "type")].loc[track_name, "parameter_value"]]
        return track_list

    def split_tracks_generator(self, resolve_attachments=True):
        if resolve_attachments:
            independent_tracks = list(self.track_metadata_df[(self.track_metadata_df["parameter"] == "type") & (self.track_metadata_df["attachment"].isna())].index)
            for track_name in independent_tracks:
                new_track = DataDF(self[["start", "end"] + list(self.track_metadata_df.loc[[track_name]]["column_name"])], parsed=True)
                new_track.track_metadata_df = self.track_metadata_df.loc[[track_name]].copy()
                new_track.track_feature_dict = {track_name: deepcopy(self.track_feature_dict[track_name])}
                new_track.track_name_list = [track_name]
                new_track.track_type_list = [self.track_metadata_df[(self.track_metadata_df["parameter"] == "type")].loc[track_name, "parameter_value"]]
                for attached_track_name in list(self.track_metadata_df[(self.track_metadata_df["parameter"] == "type") & (self.track_metadata_df["attachment"] == track_name)].index):
                    new_track.attached_track_dict[attached_track_name] = DataDF(self[["start", "end"] + list(self.track_metadata_df.loc[[attached_track_name]]["column_name"])], parsed=True)
                    new_track.attached_track_dict[attached_track_name].track_metadata_df = self.track_metadata_df.loc[[attached_track_name]].copy()
                    new_track.attached_track_dict[attached_track_name].track_feature_dict = {attached_track_name: deepcopy(self.track_feature_dict[attached_track_name])}
                    new_track.attached_track_dict[attached_track_name].track_name_list = [attached_track_name]
                    new_track.attached_track_dict[attached_track_name].track_type_list = [self.track_metadata_df[(self.track_metadata_df["parameter"] == "type")].loc[attached_track_name, "parameter_value"]]
                yield new_track
        else:
            for track_name in list(self.track_metadata_df.index.unique()):
                new_track = DataDF(self[["start", "end"] + list(self.track_metadata_df.loc[[track_name]]["column_name"])], parsed=True)
                new_track.track_metadata_df = self.track_metadata_df.loc[[track_name]].copy()
                new_track.track_feature_dict = {track_name: deepcopy(self.track_feature_dict[track_name])}
                new_track.track_name_list = [track_name]
                new_track.track_type_list = [self.track_metadata_df[(self.track_metadata_df["parameter"] == "type")].loc[track_name, "parameter_value"]]
                yield new_track

    def split_tracks_by_scaffolds(self):
        return [self.loc[scaffold_id] for scaffold_id in list(self.index.unique())]

    def split_tracks_by_scaffolds_generator(self):
        for scaffold in list(self.index.unique()):
            yield self.loc[scaffold]


if __name__ == "__main__":
    test_data = "test_data.tsv"
    test_auto_data = "test_data.auto.tsv"
    test_df = pd.read_csv(test_data, sep="\t", header=0, index_col=0, comment=None, na_values=".")
    test_auto_df = pd.read_csv(test_auto_data, sep="\t", header=0, index_col=0, comment=None, na_values=".")
    #print(test_auto_df)
    # test parsing with preset track types
    #DataDF(test_df)

    a = DataDF(test_auto_df).parse()
    print(a)

    #for da in a.split_tracks_by_scaffolds():
    #    print(da)
    #    print(da.track_metadata_df)
    #print(a.df)
    #print(a.df["admixture"]["c1"])
