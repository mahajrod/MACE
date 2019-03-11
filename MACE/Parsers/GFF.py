#!/usr/bin/env python
"""
GFF Parser Module based on pandas
"""
__author__ = 'Sergei F. Kliver'
from copy import deepcopy
from collections import OrderedDict, Iterable

import numpy as np
import pandas as pd


class CollectionGFF:

    def __init__(self, in_file=None, records=None, format="gff", parsing_mode="only_coordinates",
                 black_list=(), white_list=()):

        self.formats = ["gff", "gtf", "bed"]
        self.GFF_COLS = OrderedDict({
                                     "scaffold": 0,
                                     "source":   1,
                                     "featuretype": 2,
                                     "start": 3,
                                     "end": 4,
                                     "score": 5,
                                     "strand": 6,
                                     "phase": 7,
                                     "attributes": 8
                                     })

        self.BED_COLS = OrderedDict({
                                     "scaffold": 0,
                                     "start": 1,
                                     "end": 2
                                     })

        self.parsing_parameters = {"gff": {
                                           "only_coordinates": {
                                                                "col_names": ["scaffold", "start", "end"],
                                                                "cols":   [0, 3, 4],
                                                                "index_cols": "scaffold",
                                                                "converters": {
                                                                               "scaffold": str,
                                                                               "start":    lambda x: np.int32(x) - 1,
                                                                               "end":      np.int32,
                                                                               },
                                                                "col_name_indexes": {
                                                                                     "scaffold": 0,
                                                                                     "start":    1,
                                                                                     "end":      2
                                                                                     },
                                                                },
                                           "coordinates_and_type": {
                                                                    "col_names": ["scaffold", "featuretype", "start", "end"],
                                                                    "cols":      [0, 2, 3, 4],
                                                                    "index_cols": ["scaffold"],
                                                                    "converters": {
                                                                                   "scaffold":       str,
                                                                                   "featuretype":    str,
                                                                                   "start":          lambda x: np.int32(x) - 1,
                                                                                   "end":            np.int32,
                                                                                   },
                                                                    "col_name_indexes": {
                                                                                         "scaffold":      0,
                                                                                         "featuretype":   1,
                                                                                         "start":         2,
                                                                                         "end":           3
                                                                                         },
                                                                    },
                                           "coord_and_attr": {
                                                              "col_names": ["scaffold", "featuretype", "start", "end",
                                                                            "strand", "attributes"],
                                                              "cols":      [0, 2, 3, 4, 6, 8],
                                                              "index_cols": ["scaffold", "featuretype"],
                                                              "converters": {
                                                                             "scaffold":       str,
                                                                             "featuretype":    str,
                                                                             "start":          lambda x: np.int32(x) - 1,
                                                                             "end":            np.int32,
                                                                             "strand":         str,
                                                                             "attributes":     str,
                                                                             },
                                                              "col_name_indexes": {
                                                                                   "scaffold":      0,
                                                                                   "featuretype":   1,
                                                                                   "start":         2,
                                                                                   "end":           3,
                                                                                   "strand":        6,
                                                                                   "attributes":    8,
                                                                                   },
                                                              },
                                           "all": {
                                                   "col_names": ["scaffold", "source", "featuretype", "start", "end",
                                                                 "score", "strand", "phase", "attributes"],
                                                   "cols":      [0, 1, 2, 3, 4, 5, 6, 7, 8],
                                                   "index_cols": ["scaffold"],
                                                   "converters": {
                                                                  "scaffold":       str,
                                                                  "source":         str,
                                                                  "featuretype":    str,
                                                                  "start":          lambda x: np.int32(x) - 1,
                                                                  "end":            np.int32,
                                                                  "score":          str,
                                                                  "strand":         str,
                                                                  "phase":          str,
                                                                  "attributes":     str,
                                                                  },
                                                   "col_name_indexes": {
                                                                        "scaffold":     0,
                                                                        "source":       1,
                                                                        "featuretype":  2,
                                                                        "start":        3,
                                                                        "end":          4,
                                                                        "score":        5,
                                                                        "strand":       6,
                                                                        "phase":        7,
                                                                        "attributes":   8
                                                                        },
                                                   },
                                           "complete": {
                                                   "col_names": ["scaffold", "source", "featuretype", "start", "end",
                                                                 "score", "strand", "phase", "attributes"],
                                                   "cols":      [0, 1, 2, 3, 4, 5, 6, 7, 8],
                                                   "index_cols": ["scaffold"],
                                                   "converters": {
                                                                  "scaffold":       str,
                                                                  "source":         str,
                                                                  "featuretype":    str,
                                                                  "start":          lambda x: np.int32(x) - 1,
                                                                  "end":            np.int32,
                                                                  "score":          str,
                                                                  "strand":         str,
                                                                  "phase":          str,
                                                                  "attributes":     str,
                                                                  },
                                                   "col_name_indexes": {
                                                                        "scaffold":     0,
                                                                        "source":       1,
                                                                        "featuretype":  2,
                                                                        "start":        3,
                                                                        "end":          4,
                                                                        "score":        5,
                                                                        "strand":       6,
                                                                        "phase":        7,
                                                                        "attributes":   8
                                                                        },
                                                   },
                                           },
                                   "bed": {
                                           "only_coordinates": {
                                                                "col_names": ["scaffold", "start", "end"],
                                                                "cols":      [0, 1, 2],
                                                                "index_cols": "scaffold",
                                                                "converters": {
                                                                               "scaffold":  str,
                                                                               "start":     np.int32,
                                                                               "end":       np.int32,
                                                                               },
                                                                "col_name_indexes": {
                                                                                     "scaffold": 0,
                                                                                     "start":    1,
                                                                                     "end":      2
                                                                                     },
                                                                },
                                           }
                                   }
        self.parsing_mode = parsing_mode
        self.type_parsing_modes = ["coordinates_and_type", "all", "coord_and_attr"]
        self.attributes_parsing_modes = ["complete", "coord_and_attr"]
        self.format = format
        self.black_list = black_list
        self.white_list = white_list

        # attributes type conversion parameters
        self.parameter_separator_dict = OrderedDict()
        self.default_replace_dict = OrderedDict({
                                                 ".": None
                                                 })
        self.converters = OrderedDict()
        self.pandas_int_type_correspondence = OrderedDict({
                                                           "Int8":  np.float16,
                                                           "Int16": np.float16,
                                                           "Int32": np.float32,
                                                           "Int64": np.float64,
                                                           })

        # init aliases
        self.record_id_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["scaffold"]
        self.record_start_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["start"]
        self.record_end_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["end"]

        self.col_names = self.parsing_parameters[self.format][self.parsing_mode]["col_names"]
        self.index_cols = self.parsing_parameters[self.format][self.parsing_mode]["index_cols"]

        # load records
        if in_file:
            self.read(in_file, format=format, parsing_mode=parsing_mode, black_list=black_list, white_list=white_list)

        else:
            self.records = records

        self.scaffold_list = self.records.index.unique().to_list()

    def read(self, in_file, format="gff", parsing_mode="only_coordinates", sort=False,
             black_list=(), white_list=()):
        if format not in self.parsing_parameters:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing!" % parsing_mode)
        elif parsing_mode not in self.parsing_parameters[format]:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing in this mode(%s)!" % (format, parsing_mode))

        self.records = pd.read_csv(in_file, sep='\t', header=None, na_values=".",
                                   comment="#",
                                   usecols=self.parsing_parameters[format][parsing_mode]["cols"],
                                   converters=self.parsing_parameters[format][parsing_mode]["converters"],
                                   names=self.parsing_parameters[format][parsing_mode]["col_names"],
                                   index_col=self.parsing_parameters[format][parsing_mode]["index_cols"])

        self.records.index = pd.MultiIndex.from_arrays([self.records.index, np.arange(0, len(self.records))],
                                                       names=("scaffold", "row"))
        if parsing_mode in self.attributes_parsing_modes:
            attributes = self.parse_attributes()
            self.records.columns = pd.MultiIndex.from_arrays([
                                                              self.records.columns,
                                                              self.records.columns,
                                                              ])
            retained_columns = deepcopy(self.parsing_parameters[self.format][self.parsing_mode]["col_names"])
            for entry in "attributes", "scaffold":
                retained_columns.remove(entry)
            print self.records
            print "AAAA"
            print attributes
            self.records = pd.concat([self.records[retained_columns], attributes], axis=1)
        if sort:
            self.records.sort_values(by=["scaffold", "start", "end"])

    def parse_column(self, column, param):
        #col.replace(self.default_replace_dict, inplace=True)
        if param not in self.converters:
            return column
        elif self.converters[param] == str:
            return column
        if self.converters[param] in self.pandas_int_type_correspondence:
            col = column.apply(self.pandas_int_type_correspondence[self.converters[param]]).astype(self.converters[param])
        else:
            col = column.apply(self.converters[param])

        return col

    def parse_attributes(self):
        print("Parsing attribute field...")
        tmp_attr = self.records["attributes"].str.split(";", expand=True)
        tmp_attr_list = [tmp_attr[column].str.split("=", expand=True) for column in tmp_attr.columns]

        del tmp_attr
        attr_df_list = []

        parameter_set = set()
        for dataframe in tmp_attr_list:
            parameter_set |= set(dataframe[0].unique())

        for param in parameter_set:
            print ("\tParsing %s..." % param)
            temp_list = []
            for dataframe in tmp_attr_list:
                print("\t\tParsing fragments...")
                column = dataframe[dataframe[0] == param][1]
                if column.empty:
                    continue

                column = self.parse_column(column, param)
                temp_list.append(column)
            if not temp_list:
                continue
            print("\t\tMerging fragments...")
            tmp = pd.concat(temp_list)
            print("\t\tMerging finished...")
            del temp_list
            # TODO: check if 3 lines below are redundant in all possible cases
            shape = np.shape(tmp)
            column_number = 1 if len(shape) == 1 else shape[1]
            tmp = tmp.to_frame(param) if isinstance(tmp, pd.Series) else tmp
            tmp.columns = pd.MultiIndex.from_arrays([
                                                     ["attributes"] * column_number,
                                                     [param] * column_number
                                                     ])

            attr_df_list.append(tmp)
        #print attr_df_list
        print("Merging parameters...")
        attr = pd.concat(attr_df_list, axis=1)
        print("Merging finished.")
        attr.sort_index(level=1, inplace=True)
        print("Parsing attribute finished.")
        return attr

    def total_length(self):
        return np.sum(self.records['end'] - self.records['start'])

    def collapse_records(self, sort=True, verbose=True):
        """
        strand-independent collapse
        :param sort:
        :param verbose:
        :return:
        """
        records_before_collapse = len(self.records)

        if sort:
            self.records.sort_values(by=["scaffold", "start", "end"])
        row_list = []
        for scaffold in self.scaffold_list:
            #print scaffold
            # check if there is only one record per scaffold, necessary as pandas will return interger instead of Series
            if len(self.records.loc[[scaffold]]) == 1:
                for row in self.records.loc[[scaffold]].itertuples(index=True):
                    row_list.append(list(row))
                continue
            #print self.records.loc[scaffold]
            # remove nested records
            end_diff = self.records.loc[[scaffold]]['end'].diff()
            #print len(end_diff)
            end_diff[0] = 1
            no_nested_records_df = self.records.loc[[scaffold]][end_diff > 0]
            #print len(no_nested_records_df)
            # collapse overlapping records

            row_iterator = no_nested_records_df.itertuples(index=True)

            prev_row = list(row_iterator.next())

            for row in row_iterator:
                row_l = list(row)
                if row_l[self.record_start_col] > prev_row[self.record_end_col]:
                    row_list.append(prev_row)
                    prev_row = row_l
                else:
                    prev_row[self.record_end_col] = row_l[self.record_end_col]

            row_list.append(prev_row)
        self.records = pd.DataFrame.from_records(row_list, columns=self.col_names, index=self.index_cols)

        if verbose:
            print("Records before collapsing: %i\nRecords after collapsing: %i" % (records_before_collapse,
                                                                                   len(self.records)))

    def remove_small_records(self, min_record_length):
        records_before_collapse = len(self.records)
        self.records = self.records[(self.records['end'] - self.records['start']) >= min_record_length]
        print("Records before filtering: %i\nRecords afterfiltering: %i" % (records_before_collapse,
                                                                                   len(self.records)))

    def __add__(self, other):
        new_gff_record = CollectionGFF(records=pd.concat([self.records, other.records]),
                                       in_file=None, format=self.format,
                                       parsing_mode=self.parsing_mode,
                                       black_list=self.black_list, white_list=self.white_list)
        new_gff_record.records = new_gff_record.records.sort_values(by=["scaffold", "start", "end"])

        return new_gff_record

    def __radd__(self, other):
        new_gff_record = CollectionGFF(records=pd.concat([other.records, self.records]),
                                       in_file=None, format=other.format,
                                       parsing_mode=other.parsing_mode,
                                       black_list=other.black_list, white_list=other.white_list)
        new_gff_record.records = new_gff_record.records.sort_values(by=["scaffold", "start", "end"])

        return new_gff_record


    @staticmethod
    def get_filtered_entry_list(entry_list,
                                entry_black_list=[],
                                sort_entries=False,
                                entry_ordered_list=None,
                                entry_white_list=[]):
        white_set = set(entry_white_list)
        black_set = set(entry_black_list)
        entry_set = set(entry_list)

        if white_set:
            entry_set = entry_set & white_set
        if black_set:
            entry_set = entry_set - black_set

        filtered_entry_list = list(entry_set)
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

    def sequence_generator(self, records, sequence_collection, expression=None):
        for entry in records.itertuples():
            if expression:
                if not expression(entry):
                    continue
            yield  entry[self.record_id_col], sequence_collection[entry[self.record_id_col]][entry[self.record_start_col]:entry[self.record_end_col]]

    def extract_sequences_by_type(self, sequence_collection, record_type_black_list=[], record_type_white_list=[],
                                  return_type="collection", records_parsing_type="parse"):

        if self.parsing_mode in self.type_parsing_modes:
            if return_type == "collection":
                selected_records = self.records[self.records.index.isin(record_type_white_list, level=1) & (~self.records.index.isin(record_type_black_list, level=1))]

                from MACE.Parsers.Sequence import CollectionSequence

                extracted_records = CollectionSequence()




        else:
            pass