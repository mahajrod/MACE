#!/usr/bin/env python
"""
GFF Parser Module based on pandas
"""
__author__ = 'Sergei F. Kliver'

import os
import re

from math import sqrt
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

        # init aliases
        self.record_start_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["start"]
        self.record_end_col = self.parsing_parameters[self.format][self.parsing_mode]["col_name_indexes"]["end"]
        self.col_names = self.parsing_parameters[self.format][self.parsing_mode]["col_names"]
        self.index_cols = self.parsing_parameters[self.format][self.parsing_mode]["index_cols"]

        # load records
        if in_file:
            self.read(in_file, format=format, parsing_mode=parsing_mode, black_list=black_list, white_list=black_list)

        else:
            self.records = records

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
        for scaffold in self.records:
            # remove nested records
            end_diff = self.records.loc[scaffold]['end'].diff(axis=0)
            end_diff[0] = 1
            no_nested_records_df = self.records.loc[scaffold][end_diff > 0]

            # collapse overlapping records
            row_iterator = no_nested_records_df.itertuples(index=True)

            prev_row = row_iterator.next()

            for row in row_iterator:
                if row[self.record_start_col] > prev_row[self.record_end_col]:
                    row_list.append(prev_row)
                    prev_row = row
                else:
                    prev_row[self.record_end_col] = row[self.record_end_col]
            row_list.append(prev_row)
        self.records = pd.DataFrame(row_list, columns=self.col_names, index=self.index_cols)

        if verbose:
            print("Records before collapsing: %i\nRecords after collapsing: %i" % (records_before_collapse,
                                                                                   len(self.records)))
