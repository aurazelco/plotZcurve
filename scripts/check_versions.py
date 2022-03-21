#!/usr/bin/env python3
"""
Author: Aura Zelco

Title: plotZcurve.py

- General description:
This script imports the modules needed for plotZcurve.py and prints the version. 

- Procedure:
1. imports all necessary python modules
2. prints the version

- Usage:
This script imports the modules specified and prints the version

It is run in the command line as:
check_versions.py

- List of user-defined functions:
No user-defined functions. 

- List of imported modules:
1. numpy: to create a temporary matrix which will then be saved as dataframe 
2. pandas: to create a dataframe to then input it in R
3. rpy2 and all submodules: to install the necessary packages and to import a 
custom function from R (detailed documentation in the code)
4. math: to calculate the sqaure root for the transformation matrix


- List of known/possible bugs:
1. If the modules are not installed, the script will raise an error and exit. 
2. If other modules are added in later versions of plotZcurve.py, they have to be added also here manually. 

"""
import numpy
import pandas
import rpy2
import math
print('numpy:', numpy.__version__)
print('pandas:', pandas.__version__)
print('rpy2:', rpy2.__version__)
print('math:', math.__version__)