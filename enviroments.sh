#!/bin/bash
#Dependencies Installation
#All changes are save now, but not commited

# Print the current working directory
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

conda env create --file enviroments.yaml
