#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from BCBio import GFF
from MACE.Parsers.VCF import CollectionVCF, ReferenceGenome


def list_from_str(s):
    return s.split(",")


def figsize_from_str(s):
    try:
        x, y = map(int, s.split(','))
        return x, y
    except:
        raise argparse.ArgumentTypeError("Figsize must be x,y")

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input vcf file with mutations.")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output file with rainfall plot")
parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=300,
                    help="Dpi of figure")
parser.add_argument("-f", "--size_of_figure", action="store", dest="size_of_figure", type=figsize_from_str,
                    default=(40, 40),
                    help="Size of figure in inches. X and Y values should be separated by comma. Default: 40,40")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=list_from_str,
                    default=["svg", "eps", "pdf", "png", "jpg"],
                    help="Comma-separated list of formats (supported by matlotlib) of output figure.Default: svg,eps,pdf,png,jpg")
parser.add_argument("-l", "--suptitle", action="store", dest="suptitle",
                    help="Suptitle of figure. Default: 'Rainfall plot'")
parser.add_argument("-g", "--draw_gaps", action="store_true", dest="draw_gaps",
                    help="Draw gaps, ignored if reference genome is not set. Default: False")
parser.add_argument("-r", "--reference_genome", action="store", dest="ref_genome",
                    help="Fasta file with reference genome, required to draw gaps")
parser.add_argument("-m", "--masked_regions", action="store", dest="masked_regions",
                    help="Gff file with masked regions")
parser.add_argument("-b", "--logbase", action="store", dest="logbase", default=2, type=int,
                    help="Logbase of y axis")
args = parser.parse_args()

mutations = CollectionVCF(from_file=True, in_file=args.input, dont_parse_info_and_data=True)

if args.ref_genome:
    reference_genome = ReferenceGenome(args.ref_genome)
    reference_genome.find_gaps()
else:
    reference_genome = None

if args.masked_regions:
    masked_regions = {}
    with open(args.masked_regions) as gff_fd:
        for record in GFF.parse(gff_fd):
            masked_regions[record.id] = record
else:
    masked_regions = None

mutations.rainfall_plot(args.output_prefix, single_fig=True, dpi=args.dpi, figsize=args.size_of_figure,
                        facecolor="#D6D6D6",
                        ref_genome=reference_genome, masked_regions=masked_regions, min_gap_length=10,
                        draw_gaps=args.draw_gaps, suptitle=args.suptitle,
                        gaps_color="#777777", masked_regions_color="#aaaaaa", logbase=args.logbase,
                        extension_list=args.output_formats)
