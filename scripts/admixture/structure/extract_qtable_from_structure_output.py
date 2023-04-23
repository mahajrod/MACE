#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from MACE.Functions.General import metaopen


def convert_line_to_q_line(line):
    return line.split(":")[-1].strip()

def convert_line_to_one_space_separated_line(line):
    return " ".join(line.strip().split())

def convert_line_to_clumpp_input(line):
    line_list = list(map(lambda s: s.strip(), line.split(":")))
    number_of_preceding_columns = len(line_list[0].split())

    if number_of_preceding_columns == 4:
        return line.strip()
    elif number_of_preceding_columns == 3:
        return ":".join(line_list[:-1]) + " " + line_list[0].split()[0] + " : " + line_list[-1]
    else:
        raise ValueError("ERROR!!! Unsupported number of columns preceding  Q table!")


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input file with structure report. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file with Q table. Default: stdout")
parser.add_argument("-f", "--format", action="store", dest="format", default="q_table",
                    help="Format of output file. Allowed: 'q_table' (default), clumpp_input, 'structure'")

args = parser.parse_args()

convert_function = None
if args.format == "q_table":
    convert_function = convert_line_to_q_line
elif args.format == "clumpp_input":
    convert_function = convert_line_to_clumpp_input
elif args.format == "structure":
    convert_function = convert_line_to_one_space_separated_line
else:
    raise ValueError("ERROR!!! Unknown format {0}".format(args.format))

with metaopen(args.input, "r") as in_fd, metaopen(args.output, "w") as out_fd:
    for line in in_fd:
        if "Label (%Miss) :  Inferred clusters" in line:
            break
    for line in in_fd:
        if line == "\n":
            break
        else:
            out_fd.write(convert_function(line) + "\n")
