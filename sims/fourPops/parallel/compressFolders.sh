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
# Working directory  (first argument of command line)
# Tag that is also the name of subfolder (second argument of command line)
#wrkdir=$0
#foldertag=$1;


echo "Working directory is ${wrkdir}"
echo "Foldertag is ${foldertag}"

echo "Compressing folder ${wrkdir}/${foldertag} into ${foldertag}.tar.gz"

###
### 1. Compress folder
###
tar -zcvf ${foldertag}.tar.gz ${wrkdir}/${foldertag}
