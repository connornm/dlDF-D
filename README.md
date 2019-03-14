# dlDF-D
Connor Dolan


# Code organization

This program is for running jobs and organizing data with molecular dynamic software.
The location of directories and programs are mostly self-contained, with
use of relative paths to data, a python library, and modular executables. 

# lib.py

The main data structure can be found in the file 'lib.py', which contains a class object 'data'.
The data contained in this object is organized into two dictionary structures, 'vals' and 'meta'.
'vals' can contain the values of possible quantities of interest; lengths, energies, density matrices, ect. 
'meta' contains information about what quantities are available, as well as what systems are stored in the
data structure, and even the software being used to obtain quantities.

These dictionaries can be edited and changed freely, and the only default entries are four in 'meta'
called 'system', 'software' and 'quantity'.

Suppose that we have an instance of the data object called 'dat'.
Let's say in this object we have 3 systems; an H20 monomer, a CH4 monomer and a dimer configuration between them.
We are interested in the energy ('E') and the dipole moment ('dip'), and we use both Gaussian ('gaussian') and Orca ('orca')
to compute these quantities. The python console would give us:

	$ >>> dat.meta['system']
	$ {'H20_monomer', 'CH4_monomer', 'H20_CH4_dimer'}	
	$ >>> dat.meta['software']
	$ {'gaussian', 'orca'}
	$ >>> dat.meta['quantity']
	$ {'E', 'dip'}

The values are stored in the vals dictionary object, which is organzed such that the sets in our meta-deta
can act as keys to reference a value. The convention is data.meta[system][software][quantity]

For example, suppose we were interested in the 'CH4_monomer' energy that was obtained by Gaussian.
Using our above convention we'd type:

	$ >>> data.vals['CH4_monomer']['gaussian']['E']
	$ -44.1012

In this example our answer is given in kcal/mol, but it is ultimately up to the user on which naming units as well
as naming conventions to use. This allows for convienient comparison and organization of data between different systems
and softwares. 


# bin/

Folder of executables, mainly for file type conversion and variable retrieval. 

# data/

Folder containing all data organized by file type and optionally subcategories within each file type.

# run/

Templates for job submission based on type being submitted


