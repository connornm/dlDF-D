#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks-per-node=16
#SBATCH --output=data/slurm/CH4_dimer00304.com.out
#SBATCH --mem-per-cpu=1GB
# This is a template for running com Gaussian jobs
# Set environment variables for Gaussian
source $HOME/.g09

# Directory where the work Guassian will be done and the source directory 
wrkdir=/scratch/local/$USER/CH4_dimer00304.com
srcdir=$w

# Make a temporary work directory
mkdir -p $wrkdir
cd $wrkdir

# Copy basis functions for calculating to work 
cp ${srcdir}/data/gaussian/*gbs .
#cp ${srcdir}/data/gaussian/ParamrPW86.dat .
cp ${srcdir}/data/gaussian/ParamM05* .
# Copy Gaussian input file 
cp ${srcdir}/data/gaussian/CH4_dimer00304.com.com .

# Run Gaussian
${g09root}/g09/g09 < CH4_dimer00304.com.com > ${srcdir}/data/gaussian/CH4_dimer00304.com.log

# Remove scratch
rm -rf $wrkdir
