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

    def __init__(self, in_file=None, records=None, format="gff",
                 parsing_mode="only_coordinates"):
        self.GFF_COLS = OrderedDict({"scaffold": 0,
                                     "source":   1,
                                     "featuretype": 2,
                                     "start": 3,
                                     "end": 4,
                                     "score": 5,
                                     "strand": 6,
                                     "phase": 7,
                                     "attributes": 8})

        self.BED_COLS = OrderedDict({"scaffold": 0,
                                     "start": 1,
                                     "end": 2})
        self.only_coordinates_field_names = ["scaffold", "start", "end"]
        self.only_coordinates_converters = OrderedDict({"scaffold": str,
                                                        "start": np.int32,
                                                        "end": np.int32})

        self.coordinates_and_type_field_names = ["scaffold", "source", "start", "end"]

        self.parsing_parameters = {"gff": {
                                           "only_coordinates": {
                                                                "col_names": ["scaffold", "start", "end"],
                                                                "cols":   [0, 3, 4],
                                                                "index_cols": "scaffold",
                                                                "converters": {
                                                                               "scaffold":  str,
                                                                               "start":    np.int32,
                                                                               "end":    np.int32,
                                                                               },
                                                                },
                                           "coordinates_and_type": {
                                                                    "col_names": ["scaffold", "source", "start", "end"],
                                                                    "cols":      [0, 1, 3, 4],
                                                                    "index_colss": ["scaffold", "source"],
                                                                    "converters": {
                                                                                   "scaffold":  str,
                                                                                   "source":    str,
                                                                                   "start":     np.int32,
                                                                                   "end":       np.int32,
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
                                                                               "start":    np.int32,
                                                                               "end":      np.int32,
                                                                               },
                                                                },
                                           }
                                   }




        if in_file:
            self.read(in_file, format=format, parsing_mode=parsing_mode)

        else:
            self.records = records

    def read(self, in_file, format="gff", parsing_mode="only_coordinates"):
        if format not in self.parsing_parameters:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing!" % parsing_mode)
        elif parsing_mode not in self.parsing_parameters[format]:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing in this mode(%s)!" % (format, parsing_mode))

        self.records = pd.read_csv(in_file, sep='\t', header=None, na_values=".",
                                   usecols=self.parsing_parameters[format][parsing_mode]["cols"],
                                   converters=self.parsing_parameters[format][parsing_mode]["converters"],
                                   names=self.parsing_parameters[format][parsing_mode]["col_names"],
                                   index_col=self.parsing_parameters[format][parsing_mode]["index_cols"])
