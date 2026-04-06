<p align="center">
  <h1 align="center">MACE</h1>
</p>

<p align="center">
  <i>Tool for visualization of chromosome tracks, synteny and genomic features</i>
</p>
**MACE** is a package of scripts designed for analysis and visualization of population and comparative genomics data.
It includes a flexible pandas-based VCF parser allowing to present all the VCF file (or just a necessary columns from it) and perform the desired manipulations.
Currently, there two main direction odf the development: (1) the update of the heatmap drawing engine and (2) improvement of the pipelines related to the genome rearrangement  

# Installation
As of 06.04.2026 the latest MACE version was 1.1.38.


There are two recommended ways to install MACE.
**Option 1**: use conda to get *'relatively stable'* version of the MACE
As of 06.04.2026 the latest version MACE version was 1.1. 38
```shell
mamba install -c mahajrod routoolpa mace
```

**Option 2**: install semi manually from github to get latest version of MACE.
```shell
# install MACE and RouToolPa dependencies
mamba install 'python>=3.9' 'pandas' 'scipy' 'numpy>=1.26' 'matplotlib' 'biopython' 'bcbio-gff' 'ete3' 'statsmodels' 'pyparsing' 'xmltodict' 'venn' 'xlsxwriter'

# clone RouToolPa and MACE repositories from github
git clone https://github.com/mahajrod/RouToolPa
git clone https://github.com/mahajrod/MACE

# add RouToolPa and MACE folders to PYTHONPATH variable in your .bashrc 

```




