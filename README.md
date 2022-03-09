# Plotting a Z-curve

## Background


## Input files

### Required files

### Optional files

## Installation

### Versions of languages and packages

In the table below you can find the versions used to build this script:

| Language/Module/Library | Version |
| :------------- |:-------------:| 
| Python|  3.8.10 |
| R | 4.1.1 |
| numpy | 1.20.3 |
| pandas | 1.3.3 |
| rpy2 | 3.1.0 | 
| plot3D | | 

The versions of the python modules can be checks using check_versions.py. 

### Main scripts:

The repo can be cloned in a directory on choice as:

```shell
$ git clone https://github.com/aurazelco/BINP29_Zcurve.git
````

The scripts will then be available in the local directory, and can be run. 

### Before running:

The python script needs few modules to run properly; please install the following modules before proceeding:

The python script can be run in the command line as:
1. argparse
2. math
3. numpy
4. pandas
5. rpy2

If not present, the script will raise a ModuleNotFoundError, followed by the names of modules to be installed. 

The libraries for R should be installed the first time the software is run, if not already present. 

## Usage (v1.0.0)

Below there is a description of the command, possible to be visualized in the command line with:

```shell
$ python plotZcurve.py -h

usage: plotZcurve.py [-h] -i INPUT_GENOME -Rfunc R_SCRIPT [-o OUTPUT_PLOT] [-f OUTPUT_FORMAT [OUTPUT_FORMAT ...]]

This script reads an input genome file in a FASTA format and returns the coordinates matrix to plot a Z-curve.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT_GENOME       input genome to calculate the Z-curve
  -Rfunc R_SCRIPT       R script containing the function to plot the Z-curve
  -o OUTPUT_PLOT        optional name for saving the plot, without file extension
  -f OUTPUT_FORMAT [OUTPUT_FORMAT ...]
                        optional list of formats (separated by space): example png pdf jpeg
```

There may be a FutureWarning appearing for a pandas function, depending on the operating system. At time of release, this does not constitute a problem. 

### Examples

The sample data can be found in the corresponding folder in this repo. The genomes were retrieved  as RefSeq FASTA sequences from the NCBI database, and the links are found in the table below. 

| Species | Sample file name | Link |
| :---: |:---:| :---:|
| *E. coli* | ecoli_genome.fna | [https://www.ncbi.nlm.nih.gov/assembly/GCF_000005845.2] |
| Zika virus | zika_genome.fna | [https://www.ncbi.nlm.nih.gov/nuccore/NC_012532.1]|


## Version log

Selected updates:

```
v1.0.0		First official release Zcurve - March 2022
```