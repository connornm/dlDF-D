# dlDF-D
Connor Dolan


# Code organization
Parameter fitting for dispersionless density functionals with dispersion added in.

This program is for running jobs and organizing data with molecular dynamic software.
The location of directories and programs are mostly self-contained, with
use of relative paths to data, a python library, and modular executables. 

Heavy emphasis on consistency and standardization allows for greater flexibility and
ease of manipulating data and adapting tasks to new software. 
A keyword repeatedly used all throughout this software is the variable:

	$ NAMEKEY

The NAMEKEY variables are listed in the text file called "input", and refere to a specific systems, 
such as H20_CH4, or MY_AMAZING_MOLECULE. This way instead of having to deal with many different 
file names we can always refer to the NAMEKEY reference given in input with a file descripter after.

The other variable repeatedly used is:

	$ TYPEKEY

TYPEKEY is the variable normally referred to as a file descriptor ('com' for Gaussian), 
but can also include more specific information. 
For example, if for every NAMEKEY system you want to move one atom by
some amount, instead of making all new NAMEKEYs and running into consistency issues with
other parts of the program, it is recommended that you make a new TYPEKEY (i.e.: 'moved_atom_com')
and can now easily refer to all the altered molecules with this new type.

Currently this program is set to run with SLURM, and the submit_master and master files
are not meant to be fundamental and only used as scripts to accomplish your particular task.

The library lib.py is meant to automate all the low level file-conversion and job-running
so that high level tasks may be accomplished more flexibly and easily. 

# bin/

Folder of exectuables, mainly for file type conversion and variable retrival. 

# data/

Folder containing all data organized by file type and optionally subcategories within each file type.

# run/

Templates for job submission based on type being submitted


