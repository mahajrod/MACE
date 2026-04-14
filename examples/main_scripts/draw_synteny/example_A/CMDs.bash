#!/usr/bin/env bash

echo -e "\nExample A1. Synteny between single query and reference...\n"

draw_synteny.py -i ../../../common_files/martes_synteny/ --synteny_format psl \
                -o example_A1 --output_formats png \
                --reference MFOI --query_orderlist MZIB \
                --use_original_colors --rounded --stranded --hide_track_label

echo -e "\nExample A2. Synteny between multiple queries and reference...\n"

draw_synteny.py -i ../../../common_files/martes_synteny/ --synteny_format psl \
                -o example_A2 --output_formats png \
                --reference MFOI --query_orderlist MZIB,MMAR,MFLA \
                --use_original_colors --rounded --stranded

