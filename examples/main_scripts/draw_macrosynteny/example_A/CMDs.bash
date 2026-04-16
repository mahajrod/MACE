#!/usr/bin/env bash

echo -e "\nExample A1. Basic example. Plot is controlled by files from the input/martes_data...\n"

draw_macrosynteny.py -i input/martes_data --genome_orderlist MZIB,MMAR,MFOI,MFLA  \
                     -o example_A1 \
                     --title "Macrosynteny between martes species" \
                     --syn_file_key_column 0 --syn_file_value_column 1 \
                     --chromosome_label_fontsize 7 --genome_label_fontsize 11 \
                     --min_len_threshold 1000000 \
                     --remove_same_coords_blocks --remove_nested_blocks \
                     --output_formats svg 

echo -e "\nExample A2. Plot is mostly controlled by a config file. orderlist, invertlist, targetstrandswitchlist and querystrandswitchlis are ignored...\n"

draw_macrosynteny.py -i input/martes_data --genome_orderlist MZIB,MMAR,MFOI,MFLA  \
                     --genome_config input/example_A2.conf \
                     -o example_A2 \
                     --title "Macrosynteny between martes species" \
                     --syn_file_key_column 0 --syn_file_value_column 1 \
                     --chromosome_label_fontsize 7 --genome_label_fontsize 11 \
                     --min_len_threshold 1000000 \
                     --remove_same_coords_blocks --remove_nested_blocks \
                     --output_formats svg 

                     