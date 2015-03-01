#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from MACE.Parsers.VCF import CollectionVCF
from BCBio import GFF
from Bio import SeqIO

vcf_file = "/home/mahajrod/Genetics/MCTool/example_data/PmCDA1_3d_annotated.vcf"
masking_file = "/home/mahajrod/Genetics/MCTool/example_data/LAN210_v0.10m_masked_all_not_in_good_genes.gff"

annotation_synonym_dict = {"three_prime_UTR": "3'_UTR",
                           "five_prime_UTR": "5'_UTR",
                           "snoRNA": "ncRNA",
                           "snRNA": "ncRNA"
                           }
annotation_black_list = ["gene", "region", "ARS", "long_terminal_repeat",
                         "noncoding_exon", "intron", "repeat_region", "telomere", "gene_cassette",
                         "five_prime_UTR_intron", "LTR_retrotransposon"]

gff_file = "/home/mahajrod/Genetics/MCTool/example_data/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"


annotations_dict = SeqIO.to_dict(GFF.parse(gff_file))

annotated_mutations = CollectionVCF(from_file=True, in_file=vcf_file)
"""
annotated_mutations.get_location(annotations_dict, use_synonym=True, synonym_dict=annotation_synonym_dict)
annotated_mutations.location_pie(annotation_black_list=annotation_black_list, draw_percents_on_single_pie=False,
                                 combine_mixed=True)

annotated_mutations.rainfall_plot("rainfall")
"""



