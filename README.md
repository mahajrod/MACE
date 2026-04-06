# Description
**MACE** is a package of designed for analysis and visualization of population and comparative genomics data.
It includes a flexible pandas-based VCF parser allowing to present all the VCF file (or just a necessary columns from it) and perform the desired manipulations.
Currently, there two main direction odf the development: (1) the update of the heatmap drawing engine and (2) improvement of the pipelines related to the genome rearrangement  

# Installation
MACR is available via conda.

For installation just run command below. It should pull all the dependencies necessary for scripts itself, but you might need to install tools if they are required for scripts.

```shell
mamba install -c mahajrod routoolpa mace
```

or clone this and RouToolPa's repositories, add them to the PYTHONPATH environment variable and install all the dependencies manually.

