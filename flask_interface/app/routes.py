# imports the necessary modules for Flask to run
from flask import Flask,render_template, request, abort
# imports our custom module call 'app'
from app import app
# imports the secure_filename from the werkzeug module, to ensure secure transmission of files, since we have input files
from werkzeug.utils import secure_filename 

# same modules imported in the main python plotZcurve.py, excluding argaparse
import os
import re
import pandas as pd
import numpy as np
import math

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

# creates a path for a new directory, where the images will be temporarily stored
download_folder='app/static/images/'

# if the folder does not exists, it is created now
if not os.path.exists(download_folder):
  os.mkdir(download_folder)

# accepts only .fna files for upload
app.config['UPLOAD_EXTENSIONS'] = ['.fna']
# sets the upload path, so the files can be checked - please insert your own path
app.config['UPLOAD_PATH']= ''
# sets where the images will be downloaded and later retrieved from to be displayed
app.config['DOWNLOAD_PATH'] = download_folder
# sets the path where the R script is - please insert your own path
app.config['SCRIPT_PATH'] = ''

# loads the main page, where the files are uploaded
@app.route('/')
def main_page():
  return render_template('main_input.html')

# runs the functions contained in plotZcurve.py -> method=POST so the server knows to expect input
@app.route('/', methods=['POST'])
def plot_all():
  # initializes an empty dictionary
  file_dict={}
  # initializes a set with the allowed nucleotides
  bases = set(['a', 'c', 'g', 't'])
  # defines the transformation matrix
  tr_matrix = np.array([[1,1,-1,-1], [1,1,-1,-1], [1,-1,-1,1]])
  # multiplies the matrix for the square root of 3 divided by 4
  tr_matrix = tr_matrix*math.sqrt(3)/4
  # retrieves the R function
  Zcurve=get_Rfunc()
  # for each file uploaded:
  for uploaded_file in request.files.getlist('file'):
    # checks if the file exists; if so, returns the filename name and its path as a list
    file_check=check_files(uploaded_file)
    # filename
    filename=file_check[0].split('/')[-1].split('.')[0]
    # full path
    full_path=file_check[1]
    # checks if the file is in FASTA format
    checks_format(full_path)
    # reads the genome and returns the sequence in one string
    seq=reads_genome(full_path, bases)
    # calculates the GC content
    gc=round(GC_cont(seq),2)
    # adds the gc content to the dictionary under the filename key
    file_dict[filename] = [gc]
    # creates the matrix for plotting
    plot_matrix=creates_matrix(seq, tr_matrix)
    # puts together the full path to the plot directory where it will be saved
    out_name=os.path.join(app.config['DOWNLOAD_PATH'], filename)
    # executes the R function -> in this case, it saves the plot only as png
    Zcurve.plotZcurve(plot_matrix, out_name, 'png', filename)
    # retrieves the plot filename
    plot_name = 'images/' +  filename + '.png'
    # adds the plot filename to the dictionary under the same key as the gc content
    file_dict[filename].append([plot_name])
  # returns the html template, which will print the GC content and display the plot in the server, with an option to download the plot as png
  return render_template('print_results.html', file_dict=file_dict)


# checks if the uploaded file exists
def check_files(file):
  # retrieves the filename
  filename=secure_filename(file.filename)
  # if it is not an empty string
  if filename != '':
    # extracts the file extension
    file_ext=os.path.splitext(filename)[1]
    # if it is not in the allowed formats
    if file_ext not in app.config['UPLOAD_EXTENSIONS']:
      # server aborts
      abort(400)
    # otherwise
    else:
      # builds back the full path to the file, necessary to read it
      full_path = os.path.join(app.config['UPLOAD_PATH'], filename)
      # returns the filename and the path as list   
      return(filename, full_path)

# checks if the file is in FASTA format
def checks_format(genome):
  # opens the file
  with open(genome, 'r') as genome_file:
    # assign the first line of the file to a variable
    first_line = genome_file.readline()
    # checks if file is valid FASTA file or not
    if not first_line.startswith('>'):
      # server aborts
      abort(400)

# retrieves the genome sequence
def reads_genome(genome, bases):
  # initializes an empty string
  seq = ''
  # opens the file
  with open(genome, 'r') as genome_input:
    # reads the lines in the genome input file
    for line in genome_input:
      # if the line does not start with > (FASTA format)
      if not line.startswith('>'):
        # removes the newline and makes all letters lower case
        fragment = line.strip().lower()
        # checks if all nucleotides in the sequence are valid
        for base in fragment:            
          # if not, it raises an error and exits the script
          if base not in bases:
            # server aborts
            abort(400)
        # adds the line to the string after lowering the letters and removing the newline
        seq+=fragment
  # returns the genome sequence and the string to be used in the title       
  return(seq)

# defines a new function to count the GC content
def GC_cont(seq):
    # counts all matches of g or c in the string
    gc = len(re.findall('[gc]', seq))
    # calculates the length of the sequence, excluding Ns
    tot = len(re.findall('[acgt]', seq))
    # calculates the percentage using the whole sequence length
    perc_gc = gc * 100 / tot
    # returns the percentage
    return(perc_gc)

# calculates the coordinates matrix to be plotted
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

# imports the R function
def get_Rfunc():
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
  R_func_path =os.path.join(app.config['SCRIPT_PATH'], 'Zcurve_func.R')
  with open(R_func_path, 'r') as R_func:
    # reads the file containing R function given in the command line, and saves it in string
    string = R_func.read()
  # creates a custom module, Zcurve, which contains the R function -> now this can be used as a 
  # regular python module
  Zcurve = STAP(string, 'Zcurve')
  # returns the function as a python module now
  return Zcurve
  
