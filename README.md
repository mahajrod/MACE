<p align="center">
  <h1 align="center">MACE</h1>
</p>

<p align="center">
  <i>Tool for visualization of chromosome tracks, synteny and genomic features</i>
</p>

# Disclaimer

**MACE** is still in a pre-alpha stage and under heavy development.
The documentation has many gaps and scripts likely have many bugs or are not user-friendly yet. 
I highly encourage everyone to report any observed bugs, problems and suggestions in the _Issues_ tab.
Suggestion on new scripts and plots are very welcome too.
For specific requests or collaboration (especially projects on Mustelidae genomics), please email to _sergei.kliver@sund.ku.dk_ or _mahajrod@gmail.com_


# Installation
_As of 06.04.2026 the latest MACE version was **1.1.38**._

There are two recommended ways to install MACE.

**Option 1**: use conda to get *'relatively stable'* version of the MACE

This option makes MACE scrips available globally.

```shell
mamba install -c mahajrod routoolpa mace
```

**Option 2**: install semi manually from github to get latest version of MACE.

```shell
# Step1: install MACE and RouToolPa dependencies
mamba install 'python>=3.9' 'pandas' 'scipy' 'numpy>=1.26' 'matplotlib' 'biopython' \
              'bcbio-gff' 'ete3' 'statsmodels' 'pyparsing' 'xmltodict' 'venn' 'xlsxwriter'

# Step2: clone RouToolPa and MACE repositories from github
git clone https://github.com/mahajrod/RouToolPa
git clone https://github.com/mahajrod/MACE

# Step3: add RouToolPa and MACE folders to PYTHONPATH environment variable to your ~/.bashrc file         
```
_Option 2_ makes MACE scrips available locally from MACE/scripts folder.

# Important scripts

1. **draw_features.py**
2. **draw_variant_window_densities.py**
3. **draw_synteny.py**
4. **draw_macrosynteny.py**

# Documentation
_Wiki is under development_

# How to cite MACE

MACE is still in a pre-alpha stage and far from being published.
However, three already published articles on genomics of Mustelidae have significantly affected development of the MACE scripts.
Please, cite the **MACE repository** (https://github.com/mahajrod/MACE) and one or several of articles from the list below, depending on what scripts have you used.

1. **draw_features.py** or **draw_variant_window_densities.py**:  
Tomarovsky AA, Totikov AA, Bulyonkova TM, Perelman PL, Abramov AV, Serdyukova NA, ..., and Kliver S. 2026. Genomics of Sable (_Martes zibellina_) × Pine Marten (_Martes martes_) Hybridization. _Genome Biol Evol_ 18:evag018. https://doi.org/10.1093/gbe/evag018
2. **draw_macrosynteny.py**:  
Totikov AA, Tomarovsky AA, Perelman PL, Bulyonkova TM, Serdyukova NA, Yakupova AR, ..., and Kliver S. 2026. Comparative Genomics and Phylogenomics of the Mustelinae Lineage (Mustelidae, Carnivora). _Genome Biol Evol_ 18:evag014. https://doi.org/10.1093/gbe/evag014
3. **draw_synteny.py**:  
Kliver S, Houck ML, Perelman PL, Totikov A, Tomarovsky A, Dudchenko O, ..., and Koepfli K-P 2023. Chromosome-length genome assembly and karyotype of the endangered black-footed ferret (_Mustela nigripes_). _Journal of Heredity_:esad035. https://doi.org/10.1093/jhered/esad035


Happy drawing! 

Sergei Kliver
