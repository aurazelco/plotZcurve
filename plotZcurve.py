#!/usr/bin/env python3
"""
Author: Aura Zelco

Title: plotZcurve.py

- General description:
This script takes a genome as input, and plots the resulting Z-curve. 

- Procedure:
1. imports all necessary python modules
2. imports the R function, needed to generate the plots
3. reads the genome and stores it in a variable as a concatenated string
4. creates an empty dictionary (bases_freq), with the nucleotides (a,c,g,t) as keys and 0 as values
5. initiliazes a matrix needed for the transformation of the Z-curve (see README.md for more info)
6. initiliazes a new dictionary, to contain the values of the plot axes X, Y and Z
7. reads the genome in each position, and calculates the cumulative frequency for all bases
8. trasforms the frequencies for all bases according to the matrix, and stores the values for X, Y 
and Z in the coordinates dictionary
9. trasforms the dictionary into a pandas dataframe, which is then rendered into an R dataframe
10. plotZcurve function is used on the dataframe, and the plots are saved, either in the working 
directory or a user-defined filename, as PNG (default) or other formats. 

- Usage:
This script reads an input genome file in a FASTA format and returns a Z-curve plot. 

It is run in the command line as: python plotZcurve.py [-h] -i INPUT_GENOME 
                            -Rfunc R_SCRIPT [-o OUTPUT_PLOT] [-f OUTPUT_FORMAT [OUTPUT_FORMAT ...]]


-List of Python user-defined functions:
There are no Python user-defined functions. However, a custom R function is imported. A brief 
description is given here, but please refer to the R script for more details. 

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
I have added instructions in the README to install these packages before running the script. 

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

# defines a list of packages needed for the R function to run
packageNames = ['plot3D']
# imports a R package which is used ot check if the package in packageNames is installed
utils = rpackages.importr('utils')
# defines which CRAN mirror to check, commonly is 1
utils.chooseCRANmirror(ind=1)

# checks if the libraries are installed
packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]

# if there are libraries not previously installed:
if len(packnames_to_install) > 0:
    # it installs them
    utils.install_packages(StrVector(packnames_to_install))

# imports the library plot3D
plot3D=rpackages.importr('plot3D')

# reads the file containing R function given in the command line, and saves it in string
string = args.rscript.read()

# creates a custom module, Zcurve, which contains the R function -> now this can be used as a 
# regular python module
Zcurve = STAP(string, 'Zcurve')

#%% READ GENOME

# initializes an empty string
seq = ''
# reads the lines in the genome input file
for line in args.genome:
    # if the line does not start with > (FASTA format)
    if not line.startswith('>'):
        # add check if line does not contain AGCT!
        # adds the line to the string after lowering the letters and removing the newline
        seq+=line.strip().lower()

#%% CALCULATE COORDINATES

# initializes a dictionary where the frequencies are stored
bases_freq = {'a':0, 'g':0, 'c':0, 't':0}

# defines the transformation matrix
tr_matrix = np.array([[1,1,-1,-1], [1,1,-1,-1], [1,-1,-1,1]])
# multiplies the matrix for the square root of 3 divided by 4
tr_matrix = tr_matrix*math.sqrt(3)/4
# -> needed for the Z-curve calculations

# initializes an empty dictionary, to contain the list of values to be plotted on the 3 axes
coordinates = {'x':[], 'y':[], 'z':[]}

# for each index in the genome
for i in range(len(seq)):
    # adds to the value of the base at index i (seq[i]) 1 divided by the length of the sequence, 
    # corresponding to the frequency of one base in the whole sequence -> cumulative frequency
    bases_freq[seq[i]] += 1/(len(seq)+1)
    # creates a list with the current values present in bases_freq dictionary
    freq_values=list(bases_freq.values())
    # for each axis in the coordinates dictionary:
    for index,coord in enumerate(coordinates.keys()):
        # adds to each axis' list the sum of the transformed values using the corresponding row in the tr_matrix
        coordinates[coord].append(np.sum(freq_values*tr_matrix[index]))


#%% PLOT USING R USER-DEFINED FUNCTION

# creates a pandas dataframe out of the vertically stackes lists from the coordinates dictionary, and
# labels the 3 columns as the correspondent axes
py_df = pd.DataFrame(data=np.column_stack(list(coordinates.values())), columns=['X', 'Y', 'Z'])

# converts the pandas dataframe to a R dataframe
with localconverter(robjects.default_converter + pandas2ri.converter):
  r_coord = robjects.conversion.py2rpy(py_df)

# assigns a name to the R object
robjects.r.assign("r_coord", r_coord)

# description of parameters for Zcurve R function
'''plotZcurve function

    Parameters:
    r_coord: R dataframe
        dataframe containing the values for X, Y and Z to be plotted

    args.outfile: string
        filename to be used when saving the output files

    args.out_format: list
        list of all formats in which to save the plots

    Returns:
        Zcurve plots

'''

#executes the R function
Zcurve.plotZcurve(r_coord, args.outfile,args.out_format)



