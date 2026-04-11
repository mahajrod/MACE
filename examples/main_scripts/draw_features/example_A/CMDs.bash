#!/usr/bin/env bash

echo -e "\nExample A1. Simple command for BED as an input...\n"

draw_features.py -i input/common.records.bed -n input/common.scaffold.len -t bed -o example_A1 \
           --output_formats png

echo -e "\nExample A2. Hidden track labels...\n"

draw_features.py -i input/common.records.bed -n input/common.scaffold.len -t bed -o example_A2 \
           --hide_track_label \
           --output_formats png

echo -e "\nExample A3. Hidden track labels and rounded ends of chromosomes...\n"

draw_features.py -i input/common.records.bed -n input/common.scaffold.len -t bed -o example_A3 \
           --hide_track_label --rounded \
           --output_formats png

echo -e "\nExample A4. Hidden track labels, rounded ends of chromosomes, and renamed chromosomes...\n"

draw_features.py -i input/common.records.bed -n input/common.scaffold.len -t bed -o example_A4 \
           --hide_track_label --rounded  --scaffold_syn_file input/common.scaffold.syn \
           --syn_file_key_column 0 --syn_file_value_column 1 \
           --output_formats png

echo -e "\nExample A5. 'BED' with header as input...\n"

draw_features.py -i input/common.records_header.bed -n input/common.scaffold.len -t bed_with_header -o example_A5 \
           --hide_track_label --rounded  --scaffold_syn_file input/common.scaffold.syn \
           --syn_file_key_column 0 --syn_file_value_column 1 \
           --output_formats png

echo -e "\nExample A6. 'BED' with header as input but parsed as a table format (tsv file with header)...\n"

draw_features.py -i input/common.records_header.bed -n input/common.scaffold.len -t table -o example_A6 \
           --scaffold_column_name scaffold --start_column_name start --end_column_name end  \
           --hide_track_label --rounded  --scaffold_syn_file input/common.scaffold.syn \
           --syn_file_key_column 0 --syn_file_value_column 1 \
           --output_formats png

echo -e "\nExample A7. Syn file and centromere bed...\n"

draw_features.py -i input/common.records.bed -n input/common.scaffold.len -t bed \
           --centromere_bed input/common.scaffold.centromere.bed -o example_A7 \
           --hide_track_label --rounded  --scaffold_syn_file input/common.scaffold.syn \
           --syn_file_key_column 0 --syn_file_value_column 1 \
           --output_formats png

echo -e "\nExample A8. Syn file, centromere bed and orderlist...\n"

draw_features.py -i input/common.records.bed -n input/common.scaffold.len -t bed \
           --centromere_bed input/common.scaffold.centromere.bed -o example_A8 \
           --hide_track_label --rounded  --scaffold_syn_file input/common.scaffold.syn \
           --syn_file_key_column 0 --syn_file_value_column 1 \
           --scaffold_orderlist input/example_A8.scaffold.orderlist \
           --output_formats png

echo -e "\nExample A9. Syn file, centromere bed, orderlist and whitelist (A9)...\n"

draw_features.py -i input/common.records.bed -n input/common.scaffold.len -t bed \
           --centromere_bed input/common.scaffold.centromere.bed -o example_A9 \
           --hide_track_label --rounded  --scaffold_syn_file input/common.scaffold.syn \
           --syn_file_key_column 0 --syn_file_value_column 1 \
           --scaffold_orderlist input/example_A9.scaffold.orderlist \
           --scaffold_whitelist input/example_A9.scaffold.whitelist \
           --output_formats png

echo -e "\nExample A10. Syn file, centromere bed, orderlist and a different whitelist(A10)...\n"

draw_features.py -i input/common.records.bed -n input/common.scaffold.len -t bed \
           --centromere_bed input/common.scaffold.centromere.bed -o example_A10 \
           --hide_track_label --rounded  --scaffold_syn_file input/common.scaffold.syn \
           --syn_file_key_column 0 --syn_file_value_column 1 \
           --scaffold_orderlist input/example_A10.scaffold.orderlist \
           --scaffold_whitelist input/example_A10.scaffold.whitelist \
           --output_formats png
