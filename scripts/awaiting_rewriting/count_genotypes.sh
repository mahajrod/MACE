#!/usr/bin/env bash
# Usage:
# count_genotypes.sh <vcf_file> <0-based index of sample inv vcf file>
grep -vP "^#" $1 | cut -f $((10 + $2)) | sed s/:.*// | sort | uniq -c
