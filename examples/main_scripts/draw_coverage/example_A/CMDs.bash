#!/usr/bin/env bash

echo -e "\nExample A1. Use a table with median and mean coverages to show both...\n"

draw_coverage.py -i input/nemertina_coverage.table.gz \
                 -t table -o example_A1 \
                 -n input/nemertina_coverage.len \
                 --scaffold_whitelist input/nemertina_coverage.whitelist \
                 --scaffold_orderlist input/nemertina_coverage.orderlist \
                 -w 100000 -s 10000 --rounded --coverage_column_name median,mean -m 21,24 \
                 --output_formats png


echo -e "\nExample A2. Use the same table to show only median coverages...\n"

draw_coverage.py -i input/nemertina_coverage.table.gz \
                 -t table -o example_A2 \
                 -n input/nemertina_coverage.len \
                 --scaffold_whitelist input/nemertina_coverage.whitelist \
                 --scaffold_orderlist input/nemertina_coverage.orderlist \
                 -w 100000 -s 10000 --rounded --coverage_column_name median -m 21 \
                 --output_formats png

echo -e "\nExample A3. Use bedgraph with median coverage track as input...\n"

draw_coverage.py -i input/nemertina_coverage.median.bedgraph \
                 -t bedgraph -o example_A3 \
                 -n input/nemertina_coverage.len \
                 --scaffold_whitelist input/nemertina_coverage.whitelist \
                 --scaffold_orderlist input/nemertina_coverage.orderlist \
                 -w 100000 -s 10000 --rounded -m 21 \
                 --output_formats png

echo -e "\nExample A4. Use bedgraph with mean coverage track as input...\n"

draw_coverage.py -i input/nemertina_coverage.mean.bedgraph \
                 -t bedgraph -o example_A4 \
                 -n input/nemertina_coverage.len \
                 --scaffold_whitelist input/nemertina_coverage.whitelist \
                 --scaffold_orderlist input/nemertina_coverage.orderlist \
                 -w 100000 -s 10000 --rounded -m 24 \
                 --output_formats png

echo -e "\nExample A5. Use bedgraph with median coverage track as input and hide track labels...\n"

draw_coverage.py -i input/nemertina_coverage.median.bedgraph \
                 -t bedgraph -o example_A5 \
                 -n input/nemertina_coverage.len \
                 --scaffold_whitelist input/nemertina_coverage.whitelist \
                 --scaffold_orderlist input/nemertina_coverage.orderlist \
                 -w 100000 -s 10000 --rounded -m 21 \
                 --hide_track_label \
                 --output_formats png
