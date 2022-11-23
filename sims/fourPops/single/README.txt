# Introduction

# This folder contains some of the simulations made under a four-populations model
# It also contains all the scripts needed to run those simulations in a server with the Slurm architecture

# In this file I detail the various steps required to run these simulations

# Guidelines

# In the Slurm context, a task is to be understood as a process. 
# So a multi-process program is made of several tasks. 
# By contrast, a multithreaded program is composed of only one task, which uses several CPUs.
# Tasks are requested/created with the --ntasks option 
# while CPUs, for the multithreaded programs, are requested with the --cpus-per-task option. 

# Tasks cannot be split across several compute nodes, so requesting several CPUs 
# with the --cpus-per-task option will ensure all CPUs are allocated on the same compute node. 
# By contrast, requesting the same amount of CPUs with the --ntasks option 
# may lead to several CPUs being allocated on several, distinct compute nodes


################################################
1. login the server
################################################

# use your method of choice, normally SSH, to login the server
# You will possibly be asked to enter your passphrase.

################################################
2. GO TO THE WORKING DIRECTORY
################################################

# navigate to your working directory folder, if needed

cd /users/organization/username

# you can stay in this folder and work from here

################################################
3. CREATE FOLDER STRUCTURE
################################################

# create a folder in your working directory

mkdir fourPops

# create another folder inside this one to serve as the working directory

mkdir swe

#####################################################
4. TO COPY FROM YOUR COMPUTER TO THE SERVER USE RSYNC
#####################################################

################################################
# 4.1 copy files from my computer to the working folder in the server 
# copy R files
	
rsync -avPS /home/user/fourPops/Single/sims/*.R username@server:/users/username/fourPops

# and copy sh files
	
rsync -avPS /home/user/fourPops/Single/sims/*.sh username@server:/users/username/fourPops

################################################
5. SLURM QUEUING SYSTEM
################################################

################################################
# 5.1. Simple commands can be performed directly in your terminal, like
       - mkdir
       - cd 
       - cp
       - rm

# These can be run directly in the terminal.
# However, you cannot run in your area a program, 
# e.g. YOU CANNOT (CAN NOT) DO NOT DO THIS!!!
	
./NGSadmix

################################################
# 5.2. TO RUN a program, you need to include it in a script (*.sh)
# Scripts need to have a proper syntax, and require the following information

# All scripts need to start with the following:

  #!/bin/bash     
  #SBATCH --job-name=MyFirstSlurmJob
  #SBATCH --time=0:10:0
  #SBATCH --nodes=1
  #SBATCH --ntasks-per-node=1
  #SBATCH -p fct
  #SBATCH -q cpca73032020
  # Add the rest of your program here
  echo "Run ended"


# TO ADD PERMISSIONS TO MAKE A *.SH EXECUTABLE
chmod +x *sbatch --export=ALL,foldertag='sims' --array=1-100 launchSims.sh
.sh

# ALSO, MAKE SURE THAT THE FILES ARE IN THE CORRECT FORMAT
# particularly if they were edited or created in windows OS
dos2unix *.sh


################################################
6. TO INSTALL R PACKAGES IN CLUSTER
################################################

################################################
# 6.1. create an interactive job
# TO MAKE SMALL TESTS YOU CAN START AN INTERACTIVE JOB

srun -p xxx-q xxxx --job-name "install_R_packages" --pty bash -i

# open R
	
R

# Install the MCMC package
	
install.packages("MCMCpack", repos="http://cran.r-project.org")

# install the scrm package
	
install.packages("scrm", repos = "http://cran.r-project.org")

# there was an error detected
# package ‘scrm’ is not available (for R version 3.6.0)
# to get around this error, do the following: 

install.packages("remotes")
library(remotes)
install_version("scrm")

# quit the R session 

q(n)

# BE SURE THAT YOU EXIT THE INTERACTIVE JOB AS IT WILL COUNT TIME!!!!

Ctrl+D 

# or type 
	
exit

################################################
7. TO RUN THE SCRIPTS
################################################

################################################
# 7.1. run with 1 node and 1 task - 100 runs of 1000 sims

sbatch --export=ALL,nruns_s=1,nruns_e=100,foldertag='sims' createFolders.sh

# this will create a folder named "sims" in the directory where the sh file is located
# inside the "sims" folder, multiple folders, one for each run, will be created
# any *R files present in the original folder will also be copied to the "sims" folder

################################################
# 7.2. check the queue/partition status 

squeue -u username

################################################
# 7.3. Launch the simulations

sbatch --export=ALL,foldertag='sims' --array=1-100 launchSims.sh

################################################
# 7.4. compress the folder and download

# compress the "sims" folder 
sbatch --export=ALL,wrkdir="./swe",foldertag='sims' compressFolders.sh

# to download run this in the destinaton folder
rsync -avPS username@server:/users/username/fourPops/*.gz  /home/user/fourPops/Single/sims

################################################
7.5. if you wish, you can remove all files from the server after the download

rm -rv /users/organization/username/fourPops

################################################
7.6. check elapsed time

sacct -j  jobID --format=JobID,JobName,MaxRSS,Elapsed,CPUTime,NCPUs,NNodes

replace "jobID" by the actual number

################################################
# ARRAY JOBS 
# useful info:
# https://docs.csc.fi/computing/running/array-jobs/ 
# https://slurm.schedmd.com/job_array.html