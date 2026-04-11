#!/usr/bin/env bash

# INSTALL MAVR FIRST:
# mamba -c mahajrod mavr

WINDOW_SIZE=100000
WINDOW_STEP=10000

# get table and general stats 
get_windows_stats_mosdepth_per_base_file.py -i nemertina.coverage.bed.gz -c bed \
        -w ${WINDOW_SIZE} -s ${WINDOW_STEP} -o nemertina.coverage

# output files
# nemertina.coverage.all.stat            - whole genome metrics: length, min, max, mean, median
# nemertina.coverage.per_scaffold.stat   - per_scaffold metrics: length, min, max, mean, median
# nemertina.coverage.stat                - per window metrics: length, min, max, mean, median, uncovered, uncovered fraction

# rename and bgzip nemertina.coverage.stat as nemertina.coverage.table.gz
mv nemertina.coverage.stat nemertina.coverage.table
bgzip -f nemertina.coverage.table

# get mean coverage track in the bedgraph format
zcat nemertina_coverage.table.gz | cut -f 1,2,3,6 | tail -n +2  > nemertina_coverage.mean.bedgraph

# get median coverage track in the bedgraph format
zcat nemertina_coverage.table.gz | cut -f 1,2,3,7 | tail -n +2  > nemertina_coverage.median.bedgraph