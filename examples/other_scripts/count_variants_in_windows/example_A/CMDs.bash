#!/usr/bin/env bash

echo -e "\nExample A1. .fai as length file...\n"

count_variants_in_windows.py -i input/mustela_putorius.10x.musnig.bwa.gatk4.snp.good.masked.mean49.max2.5.min0.5.hetero.vcf.gz \
                             -n input/mustela_nigripes.v2.smithsonian.raw.fasta.fai \
                             -o example_A1

echo -e "\nExample A2. .len as length file...\n"

count_variants_in_windows.py -i input/mustela_putorius.10x.musnig.bwa.gatk4.snp.good.masked.mean49.max2.5.min0.5.hetero.vcf.gz \
                             -n input/mustela_nigripes.v2.smithsonian.raw.fasta.len \
                             -o example_A2


echo -e "\nExample A3. length file with random extension...\n"

count_variants_in_windows.py -i input/mustela_putorius.10x.musnig.bwa.gatk4.snp.good.masked.mean49.max2.5.min0.5.hetero.vcf.gz \
                             -n input/mustela_nigripes.v2.smithsonian.raw.fasta.lenttt \
                             -o example_A3