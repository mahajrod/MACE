#!/usr/bin/env bash

echo -e "\nExample_B1. Simple plot..\n"
draw_variant_window_densities.py -i input/common.data.track.bed \
    -t bedgraph -l 'SNP Densities' --colormap jet  \
    -n input/common.scaffold.len --hide_track_label --rounded \
    -o example_B1 \
    --output_formats png

echo -e "\nExample_B2. Centromeres..\n"
draw_variant_window_densities.py -i input/common.data.track.bed \
    -t bedgraph -l 'SNP Densities' --colormap jet  \
    --centromere_bed input/common.scaffold.centromere.bed \
    -n input/common.scaffold.len --hide_track_label --rounded \
    -o example_B2 \
    --output_formats png

echo -e "\nExample_B3. Centromeres + Syn file...\n"
draw_variant_window_densities.py -i input/common.data.track.bed \
    -t bedgraph -l 'SNP Densities' --colormap jet  \
    --centromere_bed input/common.scaffold.centromere.bed \
    --scaffold_syn_file input/common.scaffold.syn --syn_file_key_column 0 --syn_file_value_column 1 \
    -n input/common.scaffold.len --hide_track_label --rounded \
    -o example_B3 \
    --output_formats png

echo -e "\nExample_B4. Centromeres + Syn file + Whitelist...\n"
draw_variant_window_densities.py -i input/common.data.track.bed \
    -t bedgraph -l 'SNP Densities' --colormap jet  \
    --centromere_bed input/common.scaffold.centromere.bed \
    --scaffold_syn_file input/common.scaffold.syn --syn_file_key_column 0 --syn_file_value_column 1 \
    -n input/common.scaffold.len --scaffold_whitelist input/example_B4.scaffold.whitelist \
    --hide_track_label --rounded \
    -o example_B4 \
    --output_formats png

echo -e "\nExample_B5. Centromeres + Syn file + Whitelist + Orderlist...\n"
draw_variant_window_densities.py -i input/common.data.track.bed \
    -t bedgraph -l 'SNP Densities' --colormap jet  \
    --centromere_bed input/common.scaffold.centromere.bed \
    --scaffold_syn_file input/common.scaffold.syn --syn_file_key_column 0 --syn_file_value_column 1 \
    -n input/common.scaffold.len --scaffold_whitelist input/example_B5.scaffold.whitelist \
    --scaffold_orderlist input/example_B5.scaffold.orderlist \
    --hide_track_label --rounded \
    -o example_B5 \
    --output_formats png

echo -e "\nExample_B6. Orderlist has scaffolds absent in whitelist...\n"
draw_variant_window_densities.py -i input/common.data.track.bed \
    -t bedgraph -l 'SNP Densities' --colormap jet  \
    --centromere_bed input/common.scaffold.centromere.bed \
    --scaffold_syn_file input/common.scaffold.syn --syn_file_key_column 0 --syn_file_value_column 1 \
    -n input/common.scaffold.len --scaffold_whitelist input/example_B6.scaffold.whitelist \
    --scaffold_orderlist input/example_B6.scaffold.orderlist \
    --hide_track_label --rounded \
    -o example_B6 \
    --output_formats png

echo -e "\nExample_B7. Whitelist has scaffolds absent in the orderlist...\n"
draw_variant_window_densities.py -i input/common.data.track.bed \
    -t bedgraph -l 'SNP Densities' --colormap jet  \
    --centromere_bed input/common.scaffold.centromere.bed \
    --scaffold_syn_file input/common.scaffold.syn --syn_file_key_column 0 --syn_file_value_column 1 \
    -n input/common.scaffold.len --scaffold_whitelist input/example_B7.scaffold.whitelist \
    --scaffold_orderlist input/example_B7.scaffold.orderlist \
    --hide_track_label --rounded \
    -o example_B7 \
    --output_formats png
