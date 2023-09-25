# ZNF808 Genomics Figure Notebook

[![DOI](https://zenodo.org/badge/696224255.svg)](https://zenodo.org/badge/latestdoi/696224255)

This julia project contains data and code to generate figures for:

 *Primate-specific ZNF808 is essential for pancreatic development in humans*

 Elisa De Franco*,  Nick D L Owens*,  Hossam Montaser*,  Matthew N Wakeling,  Jonna Saarimäki-Vire,  Hazem Ibrahim, Athina Triantou,  Diego Balboa,  Richard C Caswell,  Matthew B Johnson,  Sian Ellard,  Caroline F Wright, Pancreatic Agenesis Gene Discovery Consortium,  Sarah E Flanagan,  Timo Otonkoski*,  Andrew T Hattersley*,  Michael Imbeault*

Repository contains all code to analyse the epigenetic and transcriptomic impact of ZNF808 loss during *in vitro* pancreatic differentiation. Repository contains notebook and code to regenerate all transcriptomic figures. Note repository contains source data for all gene expression and quantifications of ChIP-seq at particular loci. Repository does not contain alignment files, but code is supplied to regenerate quantifications if suitable alignment files are present.

Code and data for analysis of MER11 elements and epigentic data available in folder:
[Imbeault_lab](Imebeault_lab/)

A pre-print with a more restricted analysis is available here:
medRxiv 2021.08.23.21262262; doi: https://doi.org/10.1101/2021.08.23.21262262


## Prerequistes
Julia >= 1.8, all julia packages and their versions are specified in the included Project.toml and Manifest.toml.

Additionally Python and R paackages are used via https://github.com/JuliaPy/PyCall.jl and https://github.com/JuliaInterop/RCall.jl/ respectively:

  1. R package DESeq2 for differential expression (https://bioconductor.org/packages/release/bioc/html/DESeq2.html) ensure that R is install prior to installing `RCall` and DESeq2 is installed prior to runnning analysis.
  2. This project performs gene set enrichment via the Enrichr (https://maayanlab.cloud/Enrichr/) API using https://github.com/owensnick/Enrichr.jl. The latter package, uses Python https://requests.readthedocs.io/en/latest/ to interface with EnrichrAPI. Ensure that Python and `requests` is installed prior to installing this project.


## Installation
```bash
git clone https://github.com/owensnick/ZNF808Genomics.jl
cd ZNF808Genomics.jl
julia
```
Within julia activiate the local 
```julia
] # to enter into Pkg mode
activate .
instantiate ## for first time installation
```
To regenerate figures either use jupyter notebook within `notebooks` directory or use script as follows:
```julia
 include("notebooks/znf808_transcriptomic_figures_s4.jl")
 ```
This will generate a `figures` folder and will generate all figure panels in `svg` and `png` format.