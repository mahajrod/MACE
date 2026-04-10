#!/usr/bin/env bash

echo -e "\nExample A1. Simple command for BED as an input...\n"

draw_features.py -i input/records.bed -n input/scaffold.len -t bed -o example_A1

echo -e "\nExample A2. Hidden track labels...\n"

draw_features.py -i input/records.bed -n input/scaffold.len -t bed -o example_A2 \
           --hide_track_label

echo -e "\nExample A3. Hidden track labels and rounded ends of chromosomes...\n"

draw_features.py -i input/records.bed -n input/scaffold.len -t bed -o example_A3 \
         --hide_track_label --rounded

echo -e "\nExample A4. Hidden track labels, rounded ends of chromosomes, and renamed chromosomes...\n"

draw_features.py -i input/records.bed -n input/scaffold.len -t bed -o example_A4 \
         --hide_track_label --rounded  --scaffold_syn_file input/scaffold.syn \
         --syn_file_key_column 0 --syn_file_value_column 1

echo -e "\nExample A5. 'BED' with header as input...\n"

draw_features.py -i input/records_header.bed -n input/scaffold.len -t bed_with_header -o example_A5 \
           --hide_track_label --rounded  --scaffold_syn_file input/scaffold.syn \
           --syn_file_key_column 0 --syn_file_value_column 1

echo -e "\nExample A6. 'BED' with header as input but parsed as a table format (tsv file with header)...\n"

draw_features.py -i input/records_header.bed -n input/scaffold.len -t table -o example_A6 \
           --scaffold_column_name scaffold --start_column_name start --end_column_name end  \
           --hide_track_label --rounded  --scaffold_syn_file input/scaffold.syn \
           --syn_file_key_column 0 --syn_file_value_column 1

echo -e "\nExample A7. Syn file and centromere bed...\n"

draw_features.py -i input/records.bed -n input/scaffold.len -t bed \
         --centromere_bed input/scaffold.centromere.bed -o example_A7 \
         --hide_track_label --rounded  --scaffold_syn_file input/scaffold.syn \
         --syn_file_key_column 0 --syn_file_value_column 1

echo -e "\nExample A8. Syn file, centromere bed and orderlist...\n"

draw_features.py -i input/records.bed -n input/scaffold.len -t bed \
         --centromere_bed input/scaffold.centromere.bed -o example_A8 \
         --hide_track_label --rounded  --scaffold_syn_file input/scaffold.syn \
         --syn_file_key_column 0 --syn_file_value_column 1 \
         --scaffold_orderlist input/scaffold.orderlist

echo -e "\nExample A9. Syn file, centromere bed, orderlist and whitelist_1...\n"

draw_features.py -i input/records.bed -n input/scaffold.len -t bed \
         --centromere_bed input/scaffold.centromere.bed -o example_A9 \
         --hide_track_label --rounded  --scaffold_syn_file input/scaffold.syn \
         --syn_file_key_column 0 --syn_file_value_column 1 \
         --scaffold_orderlist input/scaffold.orderlist \
         --scaffold_whitelist input/scaffold_1.whitelist

echo -e "\nExample A10. Syn file, centromere bed, orderlist and whitelist_2...\n"

draw_features.py -i input/records.bed -n input/scaffold.len -t bed \
         --centromere_bed input/scaffold.centromere.bed -o example_A10 \
         --hide_track_label --rounded  --scaffold_syn_file input/scaffold.syn \
         --syn_file_key_column 0 --syn_file_value_column 1 \
         --scaffold_orderlist input/scaffold.orderlist \
         --scaffold_whitelist input/scaffold_2.whitelist
