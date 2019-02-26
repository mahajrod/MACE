#!/usr/bin/env python
"""
GFF Parser Module based on pandas
"""
__author__ = 'Sergei F. Kliver'

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
                                                                    "col_names": ["scaffold", "source", "start", "end"],
                                                                    "cols":      [0, 1, 3, 4],
                                                                    "index_cols": ["scaffold", "source"],
                                                                    "converters": {
                                                                                   "scaffold":  str,
                                                                                   "source":    str,
                                                                                   "start":     lambda x: np.int32(x) - 1,
                                                                                   "end":       np.int32,
                                                                                   },
                                                                    "col_name_indexes": {
                                                                                         "scaffold": 0,
                                                                                         "source":   1,
                                                                                         "start":    2,
                                                                                         "end":      3
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
        self.format = format
        self.black_list = black_list
        self.white_list = white_list

        # init aliases
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
                                   usecols=self.parsing_parameters[format][parsing_mode]["cols"],
                                   converters=self.parsing_parameters[format][parsing_mode]["converters"],
                                   names=self.parsing_parameters[format][parsing_mode]["col_names"],
                                   index_col=self.parsing_parameters[format][parsing_mode]["index_cols"])

        if sort:
            self.records.sort_values(by=["scaffold", "start", "end"])

    def total_length(self):
        return np.sum(self.records['end'] - self.records['start'])

    def collapse_records(self, sort=True, verbose=True):
        records_before_collapse = len(self.records)

        if sort:
            self.records.sort_values(by=["scaffold", "start", "end"])
        row_list = []
        for scaffold in self.scaffold_list:
            print scaffold
            # check if there is only one record per scaffold, necessary as pandas will return interger instead of Series
            if len(self.records.loc[[scaffold]]) == 1:
                for row in self.records.loc[[scaffold]].itertuples(index=True):
                    row_list.append(list(row))
                continue
            print self.records.loc[scaffold]
            # remove nested records
            end_diff = self.records.loc[[scaffold]]['end'].diff()
            print len(end_diff)
            end_diff[0] = 1
            no_nested_records_df = self.records.loc[[scaffold]][end_diff > 0]
            print len(no_nested_records_df)
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
        records_before_collapse = self.records
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
