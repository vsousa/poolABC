#!/bin/bash
#SBATCH --job-name=CreateFolders
#SBATCH --time=0:10:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p fct
#SBATCH -q cpca097952021

###
### 0. SETTINGS
###
# Number of runs (first argument of command line)
#nruns=$0;
# Name of folder (second argument of command line)
#foldertag=$1;

echo "Number of runs is $(($nruns_e-$nruns_s+1))"
echo "Foldertag is ${foldertag}"

# Path to the working directory where the sims will be run
wrkdir="/users3/fculce3c/jgcarvalho/fourPops/swe"
# Path to the folder with the R scripts to use
rdir="/users3/fculce3c/jgcarvalho/fourPops"

echo "working directory is ${wrkdir}"
echo "folder with r scripts is ${rdir}"

echo "Folder where runs subfolder will be created is ${wrkdir}/${foldertag}"
echo "R scripts will be copied from ${rdir} into ${wrkdir}/${foldertag}"

###
### 1. CREATING THE FOLDER STRUCTURE
###
# go to working directory
cd ${wrkdir}

# Create folder
mkdir ${foldertag}
cd ${foldertag}

# Create subdirectories for nruns
for (( i=${nruns_s}; i<=${nruns_e}; i++ ))
do 
	mkdir run${i}; 
done

###
### 2. COPY REQUIRED R FILES
###    copied from home folder to working directory
cp ${rdir}/*.R ${wrkdir}/${foldertag} 
