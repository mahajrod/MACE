#!/usr/bin/env bash

for SAMPLE in `bcftools query -l $1`
  do
    echo "Handling ${SAMPLE}..."
    bcftools view  --exclude-uncalled --exclude-types ref  -s ${SAMPLE} -O v $1 > $2.${SAMPLE}.vcf
  done
