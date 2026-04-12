#!/usr/bin/env bash

echo -e "\nExample A1. Basic scenario...\n"

draw_synteny.py -i ../../../common_files/martes_synteny/ --synteny_format psl \
                -o example_A1 --output_formats png \
                --reference MFOI --query_orderlist MZIB,MMAR,MFLA \
                --use_original_colors --rounded --stranded

