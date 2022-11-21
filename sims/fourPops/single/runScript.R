# read files with the necessary functions
source("functions4Pops.R")

# read the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# run with the repetition of the sims
# results are saved in folder paste("run", runi, sep="")
runi <- as.numeric(args[1]) 

# set the mean and variance
# the mean used here is the mean depth of coverage for the entire contigs of four Swedish populations - Arsklovet and Ramsö
sMean <- c(84.34, 66.76, 65.69, 68.83)
# variance for the entire contigs of four Swedish populations - Arsklovet and Ramsö
sVars <- c(1437.22, 912.43, 848.02, 1028.23) 

# create a list containing the information of the pool sizes by population
size <- rep(list(rep(5, 20)), 4)

# define how many simulations to run in this batch
nSims <- 1000

# run one batch of simulations
sims <- try(t(replicate(n = nSims, unlist(
	poolStats(model="Single", nDip=400, nPops=4, nLoci=300, nSites=2000, mutrate=1.5e-8, size=size, mean=sMean, variance=sVars, 
	          minimum=50, maximum=150, min.minor=2, NE=c(25000, 25000), ratio=c(0.1, 3), pool=c(5, 250), seq=c(0.0001, 0.001), 
	          split=c(0, 3), CW=c(1e-13, 1e-3), WC=c(1e-13, 1e-3), bT=c(0, 0.5))))))

# save the result of the simulations
print(paste("Working directory is:", getwd()))
saveRDS(sims, paste("./run", runi, "/sim", runi,".rds",sep = ""))
