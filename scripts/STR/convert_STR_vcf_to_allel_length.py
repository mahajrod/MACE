#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

import pandas as pd
import numpy as np
from RouToolPa.Parsers.VCF import CollectionVCF


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input vcf with STR calls.")
parser.add_argument("-l", "--amplicon_len_file", action="store", dest="amplicon_len_file",
                    help="Tab-separated file containing loci ids and length of amplicons calculated for reference."
                         "Optional")
parser.add_argument("--amplicon_id_column_name", action="store", dest="amplicon_id_column_name", default="primer_pair",
                    help="Column name in amplicon length file (--amplicon_len_file) to be used as amplicon ids."
                         "Ignored if --amplicon_len_file is not set.")
parser.add_argument("--amplicon_len_column_name", action="store", dest="amplicon_len_column_name",
                    default="amplicon_len",
                    help="Column name  in amplicon length file (--amplicon_len_file) to be used as amplicon length."
                         "Ignored if --amplicon_len_file is not set.")
parser.add_argument("--pop_file", action="store", dest="pop_file",
                    help="Two column tab-separated file. First column should contain sample ids, "
                         "second - population ids. Optional. Ignored if --add_population_column is not set")
parser.add_argument("--add_population_column", action="store_true", dest="add_population_column",
                    help="Add column with population ids to the output. "
                         " If --pop_file was set, population ids will be extracted from the corresponding file."
                         " Otherwise, each sample will be assigned to a separated population. ")
parser.add_argument("--encode_ids", action="store_true", dest="encode_ids",
                    help="Encode sample and population ids as numbers")
parser.add_argument("--na_replacement", action="store", dest="na_replacement", default="-100",
                    help="Value to use in output file to replace NA. Default: '-100'")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file with length of STR allels.")
parser.add_argument("-p", "--pop_df_file", action="store", dest="pop_df_file",
                    help="File to save encoded sample and pop ids. Default: not set")

args = parser.parse_args()

str_vcf_col = CollectionVCF(in_file=args.input, parsing_mode="complete")
sys.stderr.write(str(str_vcf_col))
#sys.stderr.write(str(list(str_vcf_col.records.columns)))
#sys.stderr.write("\n")
#sys.stderr.write(str(str_vcf_col.parsing_mode))
#sys.stderr.write("\n")
#sys.stderr.write(str(str_vcf_col.parsing_parameters[str_vcf_col.parsing_mode]) + "\n")
if args.pop_file:
    pop_df = pd.read_csv(args.pop_file, sep="\t", header=None, index_col="id", names=["id", "pop_id"])
else:
    pop_df = pd.DataFrame.from_records(zip(str_vcf_col.samples,
                                           range(1, len(str_vcf_col.samples) + 1)), index=["id"], columns=["id", "pop_id"])
sample_orderlist = list(pop_df.index)

pop_df["id_code"] = range(1, len(pop_df) + 1)
pop_df["pop_id_code"] = range(1, len(pop_df) + 1)

if args.amplicon_len_file:
    str_len_df = pd.read_csv(args.amplicon_len_file, sep="\t", header=0, usecols=(args.amplicon_id_column_name,
                                                                                  args.amplicon_len_column_name),
                             index_col=args.amplicon_id_column_name,)

tmp = str_vcf_col.records.reset_index(drop=False).set_index(("ID", "ID", "ID"))
if args.amplicon_len_file:
    tmp["length"] = str_len_df

genotype_df = tmp.xs("GT", level=1, axis=1)
genotype_df = genotype_df[sample_orderlist]
#genotype_df = genotype_df[sample_orderlist]
#sys.stderr.write("\n" + str(genotype_df) + "\n")
#sys.stderr.write("\n" + str(genotype_df.columns.get_level_values(0)[::-1]) + "\n")
#sys.stderr.write("\n" + str(genotype_df[genotype_df.columns.get_level_values(0)[::-1]]) + "\n")
#sys.stderr.write("\n" + str(genotype_df[genotype_df.columns.get_index_level(0)]) + "\n")

allel_length_df = tmp[["REF", "ALT"]].applymap(lambda value: len(value) if isinstance(value, str) else np.nan).astype('Int32')

if args.amplicon_len_file:
    allel_length_df = allel_length_df[["REF", "ALT"]].subtract(allel_length_df[("REF", "REF", 0)], axis=0).astype('Int32').add(tmp["length"], axis=0)

allel_length_df = allel_length_df.droplevel(0, axis=1)

max_allel_number = len(allel_length_df.columns)
genotype_column_number = len(genotype_df.columns)
merged_df = pd.merge(allel_length_df, genotype_df, left_index=True, right_index=True)
#sys.stderr.write("\n" + str(merged_df) + "\n")
genotype_length_array = []
for row in merged_df.itertuples(index=False):
    new_row = []
    for genotype in row[max_allel_number:]:
        new_row.append(row[genotype] if not pd.isnull(genotype) else np.NaN)
    genotype_length_array.append(new_row)

new_genotype_df = pd.DataFrame.from_records(genotype_length_array,
                                            index=genotype_df.index,
                                            columns=genotype_df.columns).astype('Int32')

allel_df = new_genotype_df.transpose().droplevel(level=1)
allel_df.index.name = "sample"

if args.encode_ids and args.pop_df_file:
    pop_df.to_csv(args.pop_df_file, sep="\t", header=True, index=True)

if args.encode_ids and args.add_population_column:
    allel_df["pop_id"] = pop_df["pop_id_code"]
elif args.add_population_column:
    allel_df["pop_id"] = pop_df["pop_id"]
if args.add_population_column:  # make pop_id column first
    allel_df = allel_df[[list(allel_df.columns)[-1]] + list(allel_df.columns)[:-1]]

#sys.stdout.write("\n" + str(allel_df) + "\n ")
#sys.stdout.write("\n" + str(allel_df.reindex(sample_orderlist)) + "\n ")

if args.encode_ids:
    allel_df.rename(index=pop_df["id_code"], inplace=True)
    #allel_df.sort_values(by='sample')

allel_df.to_csv(args.output, sep="\t", index=True, na_rep=args.na_replacement)
