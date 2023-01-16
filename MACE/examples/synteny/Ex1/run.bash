#!/usr/bin/env bash

SCRIPT_DIR="../../../../scripts/"
PSL_DIR="../data/psl/"
SYN_DIR="../data/syn/"
# Draw figure for psl files of 50K scale first to generate and save color scheme for chromosomes
for LEN in 50000; 
  do 
    ${SCRIPT_DIR}/draw_synteny.py -i ${PSL_DIR}/martes_martes.to.martes_foina.${LEN}.psl.gz,${PSL_DIR}/martes_zibellina.to.martes_foina.${LEN}.psl.gz \
                                  --query_labels MMAR,MZIB \
                                  --query_scaffold_white_lists ${SYN_DIR}/mmar.min_150.pseudohap2.1_HiC.purged.whitelist,${SYN_DIR}/mzib.min_150.pseudohap2.1_HiC.purged.whitelist \
                                  --reference_label MFOI \
                                  --reference_scaffold_white_list ${SYN_DIR}/mfoi.min_150.pseudohap2.1_HiC.whitelist \
                                  --query_scaffold_syn_files ${SYN_DIR}/mmar.min_150.pseudohap2.1_HiC.purged.syn,${SYN_DIR}/mzib.min_150.pseudohap2.1_HiC.purged.syn \
                                  --syn_file_key_column 0 \
                                  --syn_file_value_column 1 \
                                  -z  ${SYN_DIR}/mfoi.min_150.pseudohap2.1_HiC.orderedlist \
                                  -n ${SYN_DIR}/mfoi.min_150.pseudohap2.1_HiC.len  \
                                  --subplots_adjust_left  0.2  \
                                  --query_scaffold_order_list \
                                  ${SYN_DIR}/mmar.min_150.pseudohap2.1_HiC.purged.orderedlist,${SYN_DIR}/mzib.min_150.pseudohap2.1_HiC.purged.orderedlist  \
                                  --reference_scaffold_syn ${SYN_DIR}/mfoi.min_150.pseudohap2.1_HiC.syn  \
                                  --figure_height 0.8 \
                                  --rounded \
                                  --stranded  \
                                  -o syntheny.rounded.stranded.${LEN} \
                                  --x_tick_fontsize 15 \
                                  --title_fontsize 20 \
                                  --title "Synteny between MFOI(ref), MMAR and MZIB"\
                                  --invert_coordinates_for_target_negative_strand;
  done

# Draw figure for psl files of 100k and 1M scale

for LEN in 100000 1000000;
  do
    ${SCRIPT_DIR}/draw_synteny.py -i ${PSL_DIR}/martes_martes.to.martes_foina.${LEN}.psl.gz,${PSL_DIR}/martes_zibellina.to.martes_foina.${LEN}.psl.gz \
                                  --query_labels MMAR,MZIB \
                                  --query_color_filelist syntheny.rounded.stranded.50000.MMAR.chr_colors.tsv,syntheny.rounded.stranded.50000.MZIB.chr_colors.tsv \
                                  --query_scaffold_white_lists ${SYN_DIR}/mmar.min_150.pseudohap2.1_HiC.purged.whitelist,${SYN_DIR}/mzib.min_150.pseudohap2.1_HiC.purged.whitelist \
                                  --reference_label MFOI \
                                  --reference_scaffold_white_list ${SYN_DIR}/mfoi.min_150.pseudohap2.1_HiC.whitelist \
                                  --query_scaffold_syn_files ${SYN_DIR}/mmar.min_150.pseudohap2.1_HiC.purged.syn,${SYN_DIR}/mzib.min_150.pseudohap2.1_HiC.purged.syn \
                                  --syn_file_key_column 0 \
                                  --syn_file_value_column 1 \
                                  -z  ${SYN_DIR}/mfoi.min_150.pseudohap2.1_HiC.orderedlist \
                                  -n ${SYN_DIR}/mfoi.min_150.pseudohap2.1_HiC.len  \
                                  --subplots_adjust_left  0.2  \
                                  --query_scaffold_order_list \
                                  ${SYN_DIR}/mmar.min_150.pseudohap2.1_HiC.purged.orderedlist,${SYN_DIR}/mzib.min_150.pseudohap2.1_HiC.purged.orderedlist  \
                                  --reference_scaffold_syn ${SYN_DIR}/mfoi.min_150.pseudohap2.1_HiC.syn  \
                                  --figure_height 0.8 \
                                  --rounded \
                                  --stranded  \
                                  -o syntheny.rounded.stranded.${LEN} \
                                  --x_tick_fontsize 15 \
                                  --title_fontsize 20 \
                                  --title "Synteny between MFOI(ref), MMAR and MZIB" \
                                  --invert_coordinates_for_target_negative_strand;
  done

#draw with highlight of chromosomes with big inversions

for LEN in 100000 1000000;
  do
    ${SCRIPT_DIR}/draw_synteny.py -i ${PSL_DIR}/martes_martes.to.martes_foina.${LEN}.psl.gz,${PSL_DIR}/martes_zibellina.to.martes_foina.${LEN}.psl.gz \
                                  --query_labels MMAR,MZIB \
                                  --query_color_filelist syntheny.rounded.stranded.50000.MMAR.chr_colors.tsv,syntheny.rounded.stranded.50000.MZIB.chr_colors.tsv \
                                  --query_scaffold_white_lists ${SYN_DIR}/mmar.min_150.pseudohap2.1_HiC.purged.whitelist,${SYN_DIR}/mzib.min_150.pseudohap2.1_HiC.purged.whitelist \
                                  --reference_label MFOI \
                                  --reference_scaffold_white_list ${SYN_DIR}/mfoi.min_150.pseudohap2.1_HiC.whitelist \
                                  --query_scaffold_syn_files ${SYN_DIR}/mmar.min_150.pseudohap2.1_HiC.purged.syn,${SYN_DIR}/mzib.min_150.pseudohap2.1_HiC.purged.syn \
                                  --syn_file_key_column 0 \
                                  --syn_file_value_column 1 \
                                  -z  ${SYN_DIR}/mfoi.min_150.pseudohap2.1_HiC.orderedlist \
                                  -n ${SYN_DIR}/mfoi.min_150.pseudohap2.1_HiC.len  \
                                  --subplots_adjust_left  0.2  \
                                  --query_scaffold_order_list \
                                  ${SYN_DIR}/mmar.min_150.pseudohap2.1_HiC.purged.orderedlist,${SYN_DIR}/mzib.min_150.pseudohap2.1_HiC.purged.orderedlist  \
                                  --reference_scaffold_syn ${SYN_DIR}/mfoi.min_150.pseudohap2.1_HiC.syn  \
                                  --figure_height 0.8 \
                                  --rounded \
                                  --stranded  \
                                  -o syntheny.rounded.stranded.highlighted.${LEN} \
                                  --x_tick_fontsize 15 \
                                  --title_fontsize 20 \
                                  --title "Synteny between MFOI(ref), MMAR and MZIB" \
                                  --reference_highlight_file highlight.t \
                                  --invert_coordinates_for_target_negative_strand;
  done