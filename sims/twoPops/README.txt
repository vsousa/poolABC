# Introduction

# This folder contains some of the simulations made under a two-populations model
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
2. CREATE FOLDER STRUCTURE
################################################

# create a folder in your home space

mkdir twoPops

# inside this folder create one other folder

mkdir ark


#####################################################
3. TO COPY FROM YOUR COMPUTER TO THE SERVER USE RSYNC
#####################################################

################################################
# 3.1 copy files from my computer to the working folder in the server 
# copy R files
	
rsync -avPS /home/user/twoPops/ARK/sims/*.R username@server:/users/username/twoPops

# and copy sh files
	
rsync -avPS /home/user/twoPops/ARK/sims/*.sh username@server:/users/username/twoPops


################################################
4. TO RUN THE SCRIPTS
################################################

################################################
# 4.1. run with 1 node and 1 task - 50 runs of 5k sims each

sbatch --export=ALL,nruns_s=1,nruns_e=50,foldertag='sims' createFolders.sh

# this will create a folder named "sims" in the directory where the sh file is located
# inside the "sims" folder, multiple folders, one for each run, will be created
# any *R files present in the original folder will also be copied to the "sims" folder

# to run a second batch of simulations you should run a slightly different command:
# after finishing the first batch of simulations, run:

sbatch --export=ALL,nruns_s=51,nruns_e=100,foldertag='sims' createFolders.sh

################################################
# 4.2. check the queue/partition status 

squeue -u username

################################################
# 4.3. Launch the simulations

# for the first batch of simulations

sbatch --export=ALL,foldertag='sims' --array=1-50 launchSims.sh

# and for the second

sbatch --export=ALL,foldertag='sims' --array=51-100 launchSims.sh

################################################
# 4.4. compress the folder and download

# compress the "sims" folder 
sbatch --export=ALL,wrkdir="./ark",foldertag='sims' compressFolders.sh

# to download run this in the destinaton folder
rsync -avPS username@server:/users/username/twoPops/*.gz  /home/user/twoPops/ARK/sims

################################################
# ARRAY JOBS 
# useful info:
# https://docs.csc.fi/computing/running/array-jobs/ 
# https://slurm.schedmd.com/job_array.html
