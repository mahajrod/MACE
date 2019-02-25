#!/usr/bin/env python
"""
Sequence Parser Module
"""
__author__ = 'Sergei F. Kliver'


import re

from copy import deepcopy
from collections import OrderedDict

import numpy as np
import pandas as pd
from MACE.General import FileRoutines

from MACE.Parsers.GFF import CollectionGFF


class CollectionSequence:

    def __init__(self, in_file=None, records=None, format="fasta",
                 parsing_mode="parse", black_list=(), white_list=(),
                 masking=None, masking_file=None, masking_filetype="gff",
                 verbose=False):
        self.formats = ["fasta"]
        self.parsing_mode = parsing_mode
        self.seq_file = in_file
        self.seq_file_format = format
        self.white_list = white_list
        self.black_list = black_list

        if in_file:
            self.read(in_file, format=format, parsing_mode=parsing_mode,
                      black_list=black_list,  white_list=white_list, verbose=verbose)
        else:
            self.records = records

        if masking_file:
            print("Parsing masking...")
            self.masking = CollectionGFF(masking_file, format=masking_filetype, parsing_mode="only_coordinates",
                                         black_list=black_list, white_list=white_list)
        else:
            self.masking = masking

        self.seq_lengths = None
        self.length = 0        # None or pandas dataframe with seq_id as index
        self.scaffolds = None
        self.gaps = None          # None or pandas dataframe with seq_id as index

    @staticmethod
    def sequence_generator(sequence_file, format="fasta", black_list=(), white_list=(), verbose=False):
        if format == "fasta":
            with FileRoutines.metaopen(sequence_file, "r") as seq_fd:
                seq_id = None
                seq = ""
                for line in seq_fd:
                    if line[0] == ">":
                        if seq_id and (seq_id not in black_list):
                            if (not white_list) or (seq_id in white_list):
                                if verbose:
                                    print("Parsing %s" % seq_id)
                                yield seq_id, seq
                        seq_id = line[1:].split()[0]
                        seq = ""
                    else:
                        seq += line[:-1]
                else:
                    if seq_id and (seq_id not in black_list):
                        if (not white_list) or (seq_id in white_list):
                            if verbose:
                                print("Parsing %s" % seq_id)
                            yield seq_id, seq

    def reset_seq_generator(self):
        self.records = self.sequence_generator(self.seq_file, format=self.seq_file_format,
                                               black_list=self.black_list,  white_list=self.white_list)

    def read(self, seq_file, format="fasta", parsing_mode="generator", black_list=(), white_list=(),
             verbose=False):
        if format not in self.formats:
            raise ValueError("ERROR!!! This format(%s) was not implemented yet for parsing!" % parsing_mode)
        if parsing_mode == "generator":
            self.records = self.sequence_generator(seq_file, format=format,
                                                   black_list=black_list,  white_list=white_list,
                                                   verbose=verbose)
        elif parsing_mode == "parse":
            print("Parsing sequences...")
            self.records = OrderedDict()
            for seq_id, seq in self.sequence_generator(seq_file, format=format,
                                                       black_list=black_list,
                                                       white_list=white_list):
                self.records[seq_id] = seq

            return self.records

    def __len__(self):
        """
        :return: length of genome or None if genome was not parsed yet
        """
        return self.length

    def scaffold_number(self):
        """
        :return: number of regions/scaffolds/chromosomes in genome
        """
        return len(self.scaffolds)

    def get_stats_and_features(self, count_gaps=True, sort="True", min_gap_length=1):
        length_list = []
        gaps_list = []
        if self.parsing_mode == "generator":
            for seq_id, seq in self.records:
                length_list.append([seq_id, len(seq)])
                if count_gaps:
                    gaps_list.append(self.find_gaps_in_seq(seq, seq_id, min_gap_length=min_gap_length))

            self.reset_seq_generator()
        else:
            for seq_id in self.records:
                length_list.append([seq_id, len(self.records[seq_id])])
                if count_gaps:
                    gaps_list.append(self.find_gaps_in_seq(self.records[seq_id], seq_id, min_gap_length=min_gap_length))

        self.gaps = pd.concat(gaps_list)
        self.gaps.sort_values(by=["scaffold", "start", "end"])

        self.seq_lengths = pd.DataFrame.from_records(length_list, columns=("scaffold", "length"), index="scaffold")
        if sort:
            self.seq_lengths.sort_values(by=["length", "scaffold"])
        self.length = np.sum(self.seq_lengths["length"])
        self.scaffolds = self.seq_lengths.index.values

    @staticmethod
    def find_gaps_in_seq(sequence, seq_id=None, min_gap_length=1):
        """
        Finds gaps (N) in sequence
        :return: None
        """
        gap_reg_exp = re.compile("N+", re.IGNORECASE)
        gaps_list = []
        gaps = gap_reg_exp.finditer(sequence)  # iterator with
        for match in gaps:
            if (match.end() - match.start()) >= min_gap_length:
                gaps_list.append([match.start(), match.end()])
        gaps_list = pd.DataFrame(gaps_list, columns=("start", "end"))
        if seq_id:
            # add seq id as index
            gaps_list.index = pd.Index([seq_id] * len(gaps_list), name="scaffold")
        return gaps_list
