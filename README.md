# dlDF-D
Parameter fitting for dispersionless density functionals with dispersion added in.

This program is for fitting parameters to functionals compared to a benchmark (currently SAPT)
The location of directories and programs are mostly self-contained, with
use of relative paths to data, a python library, and modular executables. 

With proper preparation of the input and output files, the main executable program 'master' is to be submitted to SLURM by the 
submission script in the directory of it's location as follows: 

	$ sbatch submit_master

The 'master' program reads searches for a file called 'input' which contains a line separated list
of all the systems to be analyzed for a run. The scripts and programs here make repeated use of
the names listed there by always referring to it a specific SYSTEM, and when referencing the 
SYSTEM's files expects names such as 'SYSTEM.com' or 'SYSTEM.ener', ect ..

To run, you must have '.ener' files and '.com' files for each SYSTEM included in the input list,
as well as any parameters for functionals read by Gaussian. 





