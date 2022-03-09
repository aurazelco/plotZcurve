#!/usr/bin/env python3
"""
Author: Aura Zelco

Title: buscoParser.py

- General description:
This script takes a genome as input, and plots the resulting Z-curve. 

- Procedure:
1. imports all necessary python modules
2. 

- Usage:
This script reads an input genome file in a FASTA format and returns a Z-curve plot. 

It is run in the command line as: python plotZcurve.py [-h] -i INPUT_GENOME 
                            -Rfunc R_SCRIPT [-o OUTPUT_PLOT] [-f OUTPUT_FORMAT [OUTPUT_FORMAT ...]]


-List of user-defined functions:
There are no user-defined functions.

- List of imported modules:
1. argparse: a module which is used to input the different parameters.
2. math: to calculate the square root of 3
3. numpy: to create a temporary matrix which will then be saved as dataframe 
4. pandas: to create a dataframe to then input it in R
5. rpy2 and all submodules: to install the necessary packages and to import a 
custom function from R (deatiled documentation in the code)

- Possible bugs:
1. The script requires a series of python modules and R libraries. While there 
is a check for the R library to be imported, the python modules are not checked for. 
If for exmaple rpy2 is not installed, this will create an error and the script 
will exit. 
"""
#%% IMPORT MODULES

# Python modules
import argparse
import math
import numpy as np
import pandas as pd

# Modules for R
# imports the rpy2 module
import rpy2.robjects as robjects
# imports rpackages, used to check if necessary R packages are installed
import rpy2.robjects.packages as rpackages
# for installing the packages if missing
from rpy2.robjects.vectors import StrVector
# creates and export to Python a custom function made in R
from rpy2.robjects.packages import STAP
# needed for conversion from pandas dataframe to R datraframe
from rpy2.robjects import pandas2ri
# needed for conversion from pandas dataframe to R datraframe
from rpy2.robjects.conversion import localconverter

#%% ARGPARSE

# description of program, printed when -h is called in the command line
usage = 'This script reads an input genome file in a FASTA format and returns the coordinates matrix to plot a Z-curve. '

# creates an ArgumentParser object which has assigned name 'parser'
parser = argparse.ArgumentParser(description=usage)


# specification of input file - required
parser.add_argument(
    '-i',
    metavar = 'INPUT_GENOME',
    dest = 'genome',
    type=argparse.FileType('r'), # readable file
    required=True,   # needs to be inserted in the command line
    help="input genome to calculate the Z-curve" 
    )

# R script to be used to plot the Z-curve - required
parser.add_argument(
    '-Rfunc', 
    metavar = 'R_SCRIPT',
    dest = 'rscript',
    type=argparse.FileType('r'), # readable file
    required=True, # needs to be inserted in the command line
    help="R script containing the function to plot the Z-curve" 
    )

# output filename - optional
parser.add_argument(
    '-o', 
    metavar = 'OUTPUT_PLOT',
    dest = 'outfile',
    default='output',
    help="optional name for saving the plot, without file extension" 
    )

# types of plot formats to be produced - optional
parser.add_argument(
    '-f', 
    metavar = 'OUTPUT_FORMAT',
    dest = 'out_format',
    default=['png'],
    nargs='+', # there must be at least one argument if this flag is used
    help="optional list of formats (separated by space): example png pdf jpeg" 
    )

# returns result of parsing 'parser' to the class args
args = parser.parse_args()


#%% IMPORTING R USER-DEFINED FUNCTION

packageNames = ['plot3D']
utils = rpackages.importr('utils')
utils.chooseCRANmirror(ind=1)

packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]

# Running R in Python example installing packages:
if len(packnames_to_install) > 0:
    utils.install_packages(StrVector(packnames_to_install))

plot3D=rpackages.importr('plot3D')

string = args.rscript.read()

Zcurve = STAP(string, 'Zcurve')

#%% READ GENOME

seq = ''
for line in args.genome:
    if not line.startswith('>'):
        # add check if line does not contain AGCT!
        seq+=line.strip().lower()

#%% CALCULATE COORDINATES

bases_freq = {'a':0, 'g':0, 'c':0, 't':0}

tr_matrix = np.array([[1,1,-1,-1], [1,1,-1,-1], [1,-1,-1,1]])
tr_matrix = tr_matrix*math.sqrt(3)/4

coordinates = {'x':[], 'y':[], 'z':[]}

for i in range(len(seq)):
    bases_freq[seq[i]] += 1/(len(seq)+1)
    freq_values=list(bases_freq.values())
    for index,coord in enumerate(coordinates.keys()):
        coordinates[coord].append(np.sum(freq_values*tr_matrix[index]))


#%% PLOT USING R USER-DEFINED FUNCTION

py_df = pd.DataFrame(data=np.column_stack(list(coordinates.values())), columns=['X', 'Y', 'Z'])

with localconverter(robjects.default_converter + pandas2ri.converter):
  r_coord = robjects.conversion.py2rpy(py_df)


robjects.r.assign("r_coord", r_coord)


Zcurve.plotZcurve(r_coord, args.outfile,args.out_format)



