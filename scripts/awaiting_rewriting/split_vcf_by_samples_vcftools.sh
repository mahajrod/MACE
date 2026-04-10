#!/usr/bin/env bash

for SAMPLE in `bcftools query -l $1`
  do
    echo "Handling ${SAMPLE}..."
    vcf-subset --exclude-ref -c ${SAMPLE} $1 > $2.${SAMPLE}.vcf
  done
