#!/usr/bin/env bash

echo -e "\nExample_A1...\n"
draw_variant_window_densities.py -i mustela_putorius.10x.musnig.bwa.gatk4.snp.good.masked.mean49.max2.5.min0.5.hetero.vcf.gz  \
                                 -o example_A1 \
                                 -w 1000000 \
                                 -s 100000 \
                                 -a mustela_nigripes.v2.smithsonian.chr.whitelist  \
                                 -z mustela_nigripes.v2.smithsonian.chr.orderlist \
                                 -n mustela_nigripes.v2.smithsonian.raw.fasta.len \
                                 --scaffold_syn_file mustela_nigripes.v2.smithsonian.chr.syn  \
                                 --syn_file_key_column 1 \
                                 --syn_file_value_column 0 \
                                 --density_thresholds 0,0.1,0.5,0.75,1.0,1.25,1.5,2.0,2.5 \
                                 --hide_track_label \
                                 --rounded
