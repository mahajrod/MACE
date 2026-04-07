#!/usr/bin/env bash

echo -e "\nExample 1. Simple command for BED as an input...\n"

draw_features.py -i records.bed -n scaffold.len -t bed -o example_1

echo -e "\nExample 2. Hidden track labels...\n"

draw_features.py -i records.bed -n scaffold.len -t bed -o example_2 \
           --hide_track_label

echo -e "\nExample 3. Hidden track labels and rounded ends of chromosomes...\n"

draw_features.py -i records.bed -n scaffold.len -t bed -o example_3 \
         --hide_track_label --rounded

echo -e "\nExample 4. Hidden track labels, rounded ends of chromosomes, and renamed chromosomes...\n"

draw_features.py -i records.bed -n scaffold.len -t bed -o example_4 \
         --hide_track_label --rounded  --scaffold_syn_file scaffold.syn \
         --syn_file_key_column 0 --syn_file_value_column 1

echo -e "\nExample 5. 'BED' with header as input...\n"

draw_features.py -i records_header.bed -n scaffold.len -t bed_with_header -o example_5 \
           --hide_track_label --rounded  --scaffold_syn_file scaffold.syn \
           --syn_file_key_column 0 --syn_file_value_column 1

echo -e "\nExample 6. 'BED' with header as input but parsed as a table format (tsv file with header)...\n"

draw_features.py -i records_header.bed -n scaffold.len -t table -o example_6 \
           --scaffold_column_name scaffold --start_column_name start --end_column_name end  \
           --hide_track_label --rounded  --scaffold_syn_file scaffold.syn \
           --syn_file_key_column 0 --syn_file_value_column 1
