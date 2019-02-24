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
        self.GFF_COLUMNS = OrderedDict({"scaffold": 0,
                                        "source":   1,
                                        "featuretype": 2,
                                        "start": 3,
                                        "end": 4,
                                        "score": 5,
                                        "strand": 6,
                                        "phase": 7,
                                        "attributes": 8})

        self.BED_COLUMNS = OrderedDict({"scaffold": 0,
                                        "start": 1,
                                        "end": 2})
        self.only_coordinates_field_names = ["scaffold", "start", "end"]
        self.only_coordinates_converters = OrderedDict({"scaffold": str,
                                                        "start": np.int32,
                                                        "end": np.int32})

        self.coordinates_and_type_field_names = ["scaffold", "source", "start", "end"]

        if in_file:
            self.read(in_file, format=format, parsing_mode=parsing_mode)

        else:
            self.records = records

    def read(self, in_file, format="gff", parsing_mode="only_coordinates"):

        if parsing_mode == 'only_coordinates':
            if format == "gff":
                columns_to_read = [self.GFF_COLUMNS[entry] for entry in self.only_coordinates_field_names]
            elif format == "bed":
                columns_to_read = [self.BED_COLUMNS[entry] for entry in self.only_coordinates_field_names]
            else:
                raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing in this mode(%s)!" % (format, parsing_type))
            self.records = pd.read_csv(in_file, sep='\t', header=None, na_values=".",
                                       usecols=columns_to_read,
                                       converters=self.only_coordinates_converters,
                                       names=self.only_coordinates_field_names, index_col=self.GFF_COLUMNS["scaffold"])

        elif parsing_mode == 'coordinates_and_type':
            if format == "gff":
                columns_to_read = [self.GFF_COLUMNS[entry] for entry in self.coordinates_and_type_field_names]
            else:
                raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing in this mode(%s)!" % (format, parsing_type))
            self.records = pd.read_csv(in_file, sep='\t', header=None, na_values=".",
                                       usecols=columns_to_read,
                                       converters=self.only_coordinates_converters,
                                       names=self.coordinates_and_type_field_names,
                                       index_col=(self.GFF_COLUMNS["scaffold"],
                                                  self.GFF_COLUMNS["source"]))
