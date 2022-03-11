#!/usr/bin/env python3
"""
Author: Aura Zelco

Title: plotZcurve.py

- General description:
This script takes a genome as input, and plots the resulting Z-curve. 

- Procedure:
1. imports all necessary python modules
2. imports the R function, needed to generate the plots
3. for each genome, reads the genome and stores it in a variable as a concatenated string (if the file is indeed in FASTA format)
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

It is run in the command line as:
    plotZcurve.py [-h] -i INPUT_GENOME [INPUT_GENOME ...] [-f OUTPUT_FORMAT [OUTPUT_FORMAT ...]]
                     [-o OUTPUT_PATH] [-s SCRIPT_PATH] [-gc] [-out_gc OUTPUT_GC]


- List of user-defined functions:
1. dir_path: checkes if the directory exists
2. checks_input: checks if the genome is in FASTA format
3. reads_genome: cretaes one string from the genome sequence and extract the filename, used later
4. GC_cont: calculates the GC content in the sequence
5. creates_matrix: from the genome sequence string, creates the matrix with coordinates to be plotted 

Zcurve: a custom R function is imported; a brief description is given here, but please refer to the R script for more details. 


- List of imported modules:
1. argparse: a module which is used to input the different parameters
2. os: to retrieve the current working directory
3. re: to work with Regex
4. math: to calculate the square root of 3
5. numpy: to create a temporary matrix which will then be saved as dataframe 
6. pandas: to create a dataframe to then input it in R
7. rpy2 and all submodules: to install the necessary packages and to import a 
custom function from R (detailed documentation in the code)

- Possible errors addressed in the script:
1. InvalidInput: if the input file does not start either with > (fasta format)
2. InvalidNucleotide: if there are non-nucleotides characters in the sequence


- List of known/possible bugs:
1. At the time of release, a function from the rpy2 module can raise a FutureWarning in 
certain operating systems. However, the code still runs, so it is not addressed at the moment. 
2. Is the user inputs the -out_gc but not -gc flag, the script will still print 
the GC content on the terminal and not in the output GC filename specified
3. The script requires a series of python modules and R libraries. While there 
is a check for the R library to be imported, the python modules are not checked for. 
If for example rpy2 is not installed, this will create an error and the script 
will exit. 
I have added instructions in the README to install these packages before running the script. 
4. The input file has to be in FASTA format, and the script will not recognize if multiple genes/genomes 
are concatenated in the same file; they will be treated as one sequence. 

"""
#%% IMPORT MODULES

# Python modules
import argparse
import re
import os 
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
    required=True, 
    nargs='+', # there must be at least one argument if this flag is used
    help="input genome(s) to calculate the Z-curve, can be more than one" 
    )

# types of plot formats to be produced - optional
parser.add_argument(
    '-f', 
    metavar = 'OUTPUT_FORMAT',
    dest = 'out_format',
    default=['png'],
    nargs='+', # there must be at least one argument if this flag is used
    help="optional: list of formats (separated by space): example png pdf jpeg" 
    )

# output path - where to save output plots - optional
parser.add_argument(
    '-o', 
    metavar = 'OUTPUT_PATH',
    dest = 'out_path',
    type=os.path.abspath, # extracts the absolute path, easier to navigate through the tree
    default = os.path.curdir, # the default is the present working directory
    help="optional: path to output directory" 
    )

# R script path - in case R script is not in current working directory - optional
parser.add_argument(
    '-s', 
    metavar = 'SCRIPT_PATH',
    dest = 'script_path',
    type=os.path.abspath, # extracts the absolute path, easier to navigate through the tree
    default = os.path.curdir, # the default is the present working directory
    help="path to Zcurve_func.R, needed if the R script is not in the current working directory" 
    )

# GC content - if the user wants the GC content saved in a file instead of printed on the screen - optional
parser.add_argument(
    '-gc', 
    dest = 'save_gc',
    action="store_true",
    help="optional: in case -save_gc is used, the script will save the GC content calculations to a file instead of printing to the console" 
    )

# GC content output name - if the user wants to specify the output filename for the GC content - optional
parser.add_argument(
    '-out_gc',
    metavar = 'OUTPUT_GC',
    dest = 'out_gc',
    type=argparse.FileType('w'), 
    default='GC_content_output.txt',
    help= "optional: output file where the GC content will be written in the -gc flag is used (default 'GC_content_output.txt' in the working directory)" 
    )

# returns result of parsing 'parser' to the class args
args = parser.parse_args()


#%% CUSTOM ERRORS

'Creates a new class of custom error messages in this script'
class CustomError(Exception):
    pass

'Raised if input file is not a fasta file'
class InvalidInput(CustomError):
    pass

'Raised if a sequence contains non-nucleotide characters'
class InvalidNucleotide(CustomError):
    pass


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

R_func_path = args.script_path + '/Zcurve_func.R'

with open(R_func_path, 'r') as R_func:
    # reads the file containing R function given in the command line, and saves it in string
    string = R_func.read()

# creates a custom module, Zcurve, which contains the R function -> now this can be used as a 
# regular python module
Zcurve = STAP(string, 'Zcurve')

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

#%% USER-DEFINED FUNCTIONS - CHECKS IF OPTIONAL OUTPUT DIRECTORY EXISTS

'''DIR_PATH

    Parameters
    ----------
    seq : string
        folder path as string

    Returns
    -------
    string: folder path as string

'''

# checks if directory is valid
def dir_path(string):
    # if the directory exists, return the direcotry path as it is
    if os.path.isdir(string):
        return string
    # if not, raise an error
    else:
        raise argparse.ArgumentTypeError("{path_dir} is not a valid path")

#%% USER-DEFINED FUNCTIONS - READ GENOME

def checks_input(genome):
    # assign the first line of the file to a variable
    first_line = genome.readline()
    # checks if file is valid FASTA file or not
    if not first_line.startswith('>'):
        raise InvalidInput('Your input file {} is not valid. Please insert a fasta file' .format(genome))

def reads_genome(genome):
    # initializes an empty string
    seq = ''
    # initializes a set to check if the sequence contains other characters than nucleotides
    bases = set(['a', 'c', 'g', 't'])

    # extracts the filename to be used as title of the plot: splits by /, and retrieves the last element
    # which is going to be the name, and keeps only the name and not the file format eg '.fna'; 
    # will also be sued for the output
    plot_main=genome.name.split('/')[-1].split('.')[0]

    # reads the lines in the genome input file
    for line in genome:
        # if the line does not start with > (FASTA format)
        if not line.startswith('>'):
            # add check if line does not contain AGCT!
            fragment = line.strip().lower()
            # checks if all nucleotides in the sequence are valid
            for base in fragment:
                # if not, it raises an error and exits the script
                if base not in bases:
                    raise InvalidNucleotide('Your input file contain unvalid nucleotides. Please insert a valid input fasta file')
            # adds the line to the string after lowering the letters and removing the newline
            seq+=fragment
    # returns the genome sequence and the string to be used in the title       
    return(seq, plot_main)

''' GC_CONT

    Parameters
    ----------
    seq : string
        nucleotide sequence

    Returns
    -------
    gc_cont: float
        GC percentage in the sequence
'''

# defines a new function to count the GC content
def GC_cont(seq):
    # counts all matches of G or C in the string
    gc = len(re.findall('[gc]', seq))
    # calculates the length of the sequence, excluding Ns
    tot = len(re.findall('[acgt]', seq))
    # calculates the percentage using the whole sequence length
    perc_gc = gc * 100 / tot
    # returns the percentage
    return(perc_gc)

#%% USER-DEFINED FUNCTIONS - CALCULATE COORDINATES MATRIX

def creates_matrix(seq, tr_matrix):
    # initializes a dictionary where the frequencies are stored
    bases_freq = {'a':0, 'g':0, 'c':0, 't':0}

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

    # creates a pandas dataframe out of the vertically stackes lists from the coordinates dictionary, and
    # labels the 3 columns as the correspondent axes
    py_df = pd.DataFrame(data=np.column_stack(list(coordinates.values())), columns=['X', 'Y', 'Z'])
    # converts the pandas dataframe to a R dataframe
    with localconverter(robjects.default_converter + pandas2ri.converter):
      r_coord = robjects.conversion.py2rpy(py_df)
    # assigns a name to the R object
    robjects.r.assign("r_coord", r_coord)
    
    # returns the coordinates dictionary
    return(r_coord)


#%% MAIN

out_path=dir_path(args.out_path)

# defines the transformation matrix
tr_matrix = np.array([[1,1,-1,-1], [1,1,-1,-1], [1,-1,-1,1]])
# multiplies the matrix for the square root of 3 divided by 4
tr_matrix = tr_matrix*math.sqrt(3)/4
# -> needed for the Z-curve calculations

if args.out_gc and not args.save_gc:
    print('The GC content will be printed to the terminal. If you want to save the GC content in an output file, please add the -gc flag to the command')

# for each genome in the list provided after the -i flag
for genome_input in args.genome:
    # checks if the input is in FASTA format
    checks_input(genome_input)
    # extracts the sequence and the genome filename and saves them in a list
    params=reads_genome(genome_input)
    # assigns the first element of the list (the whole genome sequence) to seq
    seq=params[0]
    # assigns the genome filename
    file_name=params[1]
    # after the file has been read, it calculates the GC content on the whole genome
    gc_file = GC_cont(seq)
    # if the -gc flag is used
    if args.save_gc:
        # prints the filename and the GC content to the out_gc file
        print('{}: {:.2f}%' .format(file_name, gc_file), file=args.out_gc)
    # if not, prints to the terminal
    else:
        # prints the filename and the GC content to the console
        print('{}: {:.2f}%' .format(file_name, gc_file))
    # combines the output plot name 
    out_name=f'{out_path}/{file_name}'
    # creates the matrix needed to run the plotting function
    plot_matrix=creates_matrix(seq, tr_matrix)
    # executes the R function and generates the plot
    Zcurve.plotZcurve(plot_matrix, out_name, args.out_format, file_name)


