# check the installed packages and add if missing
list.of.packages <- c("MCMCpack", "scrm") # list the required packages here
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
# install missing packages
if(length(new.packages) > 0) {
  install.packages(new.packages, lib = "/users3/fculce3c/jgcarvalho/R/x86_64-redhat-linux-gnu-library/3.6", 
                   repos = "http://cran.r-project.org")
}

# load the packages
for(i in 1:length(list.of.packages)) {
  library(list.of.packages[i], character.only = TRUE)  
}

# Create parameters vector --------------------------------------------------------------

# CREATEPARAMS
createParams <- function(NE, ratio, split, pool, seq, CW = NA, WC = NA, Mb = NA, Md = NA, site = NA, anc = NA, bT = NA, 
                         bCW = NA, bWC = NA, pMd = NA, model, digits = 5) {
  
  # check if the input is correct - the selected model should be one of the following
  if(model %in% c("Single", "Parallel", "star", "noMig", "2pops", "simple2pops", "mig2pops") == FALSE) 
    stop(paste("The selected model should be either Single, Parallel, star, noMig, 2pops, simple2pops or mig2pops. Please check"))
 
  # Draw a value for the population size of the most ancestral population - it's also the Ne
  Ne <- runif(n = 1, min = NE[1], max = NE[2])
  # Draw the values for the extant and ancestral population sizes
  # Values are drawn as a ratio of the Ne size - a value of 2 means that the population is twice the size of the Ne 
  N1 <- round(runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  N2 <- round(runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  # set the values in the natural scale
  N1 <- 10^N1; N2 <- 10^N2
  
  # Draw the split times from a uniform distribution - times are drawn on a 4Ne scale 
  Split <- runif(n = 1, min = split[1], max = split[2])
  # Draw the errors for the pooling parameter
  Pool_Error <- runif(n = 1, min = pool[1], max = pool[2])
  # and the sequencing error
  Error <- runif(n = 1, min = seq[1], max = seq[2])
  
  # stop the function if we are working with the noMig model - only requires 6 parameters
  if(model == "noMig") {
    # Create the parameters vector for this particular model
    parameters <- c(Ne, N1, N2, Split, Pool_Error, Error)
    # add names to the columns
    names(parameters) <- c("Ne", "N1", "N2", "Split", "Pool_Error", "Error")
    # stop the function and output the parameters vector
    stop(return(parameters))
  }
  
  # stop the function if we are working with the mig2pops model
  if(model == "mig2pops") {
    
    # draw the baseline migration rate (at the m scale) - this is the migration between different ecotypes in the same site
    Mb <- runif(n = 1, min = Mb[1], max = Mb[2])
    # draw the migration rate that differs from the baseline migration rate - also between different ecotypes at the same site 
    Md <- runif(n = 1, min = Md[1], max = Md[2])
    
    # proportion of loci with a migration rate that differs from the baseline migration rate 
    pMd <- runif(n = 1, pMd[1], pMd[2])
    # assume that the proportion of loci with a baseline migration rate is 1 - pMd
    pMb <- 1 - pMd
    
    # Create the parameters vector for this particular model
    parameters <- c(Ne, N1, N2, Split, Pool_Error, Error, Mb, Md, pMb, pMd)
    # add names to the columns
    names(parameters) <- c("Ne", "N1", "N2", "Split", "Pool_Error", "Error", "Mb", "Md", "pMb", "pMd")
    # stop the function and output the parameters vector
    stop(return(parameters))
  }
  
  # Draw the migration rate (at the mCW scale) - 
  # this corresponds to the migration between different ecotypes in the same site
  # The migration rate from crab to wave is drawn at the mCW scale
  mCW <- runif(n = 1, min = CW[1], max = CW[2])
  # and the migration rate from wave to crab is drawn as a ratio of the migration rate from crab to wave
  mWC <- runif(n = 1, min = WC[1], max = WC[2])
  
  # proportion of the genome without migration - total barrier
  total <- rbeta(n = 1, shape1 = 1, shape2 = 10)
  # replace values below the minimum threshold with the minimum
  total[total < 1e-2] <- bT[1]
  # and values above the maximum threshold with the maximum
  total[total > bT[2]] <- bT[2]
  
  # assume that the proportion of the genome with unrestricted migration is 1 - proportion without migration
  pMig <- 1 - total
  
  # stop the function if we are working with the simple 2 pops model - only requires 10 parameters
  if(model == "simple2pops") {
    # Create the parameters vector for this particular model
    parameters <- c(Ne, N1, N2, Split, Pool_Error, Error, mCW, mWC, pMig, total)
    # add names to the columns
    names(parameters) <- c("Ne", "N1", "N2", "Split", "Pool_Error", "Error", "mCW", "mWC", "pM", "pNO")
    # stop the function and output the parameters vector
    stop(return(parameters))
  }
  
  # stop the function if we are working with the 2pops model - only requires 12 parameters
  if(model == "2pops") {
    # Create the parameters vector for this particular model
    parameters <- c(Ne, N1, N2, Split, Pool_Error, Error, mCW, mWC, pMig, pCW, pWC, total)
    # add names to the columns
    names(parameters) <- c("Ne", "N1", "N2", "Split", "Pool_Error", "Error", "mCW", "mWC", "pM", "pCW", "pWC", "pNO")
    # stop the function and output the parameters vector
    stop(return(parameters))
  }
  
  # Draw the values for the other extant populations and ancestral population sizes
  # Values are drawn as a ratio of the Ne size - a value of 2 means that the population is twice the size of the Ne 
  N3 <- round(runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  N4 <- round(runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  # set the values in the natural scale
  N3 <- 10^N3; N4 <- 10^N4
  
  # Draw additional migration rates - this corresponds to the migration between the same ecotype in different sites
  # it is also drawn as a ratio of the migration rate from crab to wave
  #mCC <- runif(n = 1, min = site[1], max = site[2])
  #mWW <- runif(n = 1, min = site[1], max = site[2])
  
  # stop the function if we are working with the star model - only requires 18 parameters
  if(model == "star") {
    # Create the parameters vector for this particular model
    parameters <- c(Ne, N1, N2, N3, N4, Split, Pool_Error, Error, mCW, mWC, mCC, mWW, pMig, pCW, pWC, total)
    # add names to the columns
    names(parameters) <- c("Ne", "N1", "N2", "N3", "N4", "Split", "Pool_Error", "Error", "mCW", "mWC", "mCC", 
                           "mWW", "pM", "pCW", "pWC", "pNO")
    # stop the function and output the parameters vector
    stop(return(parameters))
  }
  
  # Draw the values for the ancestral population sizes
  # Values are drawn as a ratio of the Ne size - a value of 2 means that the population is twice the size of the Ne 
  NA1 <- round(runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  NA2 <- round(runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  # set the values in the natural scale
  NA1 <- 10^NA1; NA2 <- 10^NA2
  
  # Draw the second split time from a uniform distribution - times are drawn on a 4Ne scale 
  Dsplit <- runif(n = 1, min = split[1], max = split[2])
  
  # and finally, the migration rates between the ancestral populations 
  # these are also drawn as a ratio of the migration rate from crab to wave
  # mAA <- runif(n = 1, min = anc[1], max = anc[2])
  
  # Create the parameters vector for the single and parallel models
  parameters <- c(Ne, N1, N2, N3, N4, NA1, NA2, Split, Dsplit, Pool_Error, Error, mCW, mWC, pMig, total)
  # add names to the columns
  names(parameters) <- c("Ne", "N1", "N2", "N3", "N4", "NA1", "NA2", "Split", "Dsplit", "Pool_Error", "Error", 
                         "mCW", "mWC", "pM", "pNO")
  # output the parameters vector
  parameters
}


# SCRM command line -------------------------------------------------------

# CMD SINGLE
cmdSingle <- function(parameters, nSites, nLoci, nDip, nPops, mutrate) {
  
  # Read the vector with the parameters and assign each parameter to the correct command name
  Ne <- parameters[1]
  # set the relative size of each population - for the extant populations
  N1 <- parameters[2]; N2 <- parameters[3]; N3 <- parameters[4]; N4 <- parameters[5]
  # and the ancient populations
  NA1 <- parameters[6]; NA2 <- parameters[7]
  
  # get the proportion of loci with migration - no barriers against migration between the different ecotypes
  pM <- parameters[14]
  # and the proportion of loci without any migration - total barrier against migration between the different ecotypes
  pNO <- parameters[15]
  
  # get the migration rates - between different ecotypes at the same site - this is the m value on the M = 4N0m formula
  # from the crab ecotype to the wave ecotype
  mCW <-  parameters[12]
  # from the wave ecotype to the crab ecotype
  mWC <-  parameters[13]
  
  # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time) 
  # and REMEMBER that M = 4N0m
  # set the migration rates to the scale of Nref - between different ecotypes at the same site
  # from the crab ecotype to the wave ecotype
  mig_CW <- 4*Ne*mCW # -m 2 1 mig_CW
  # from the wave ecotype to the crab ecotype
  mig_WC <- 4*Ne*mWC # -m 1 2 mig_WC
  
  # get the time of the recent split event
  Rsplit <- round(parameters[8], digits = 3)
  # and the ancient split event
  Asplit <- round(parameters[9], digits = 3)
  # the actual time of the ancient split event is obtained by doing:
  Asplit <- Rsplit + Asplit
  
  # Compute the value of theta
  mutrate_locus <- nSites*mutrate
  theta <- 4*Ne*mutrate_locus
  # Create a vector with information about how many haplotypes are sampled from each population
  n <- c(rep(nDip/nPops, times = nPops))
  
  # create a variable to be used for a minor correction related to the split time
  # this will ensure that the changes in migration rates occur at different times
  # correctionTime <- 1e-10
  
  # use a multinomial distribution to get the number of loci simulated under each category
  lociTotal <- as.vector(rmultinom(n = 1, size = nLoci, c(pM, pNO)))
  
  ## cheat code: pop1 - crab in site 1; pop2 - wave in site 1; pop3 - crab in site 2; pop4 - wave in site 2
  # create command line with no barriers to migration
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and 
  # nrep: the number of independent loci that will be produced
  with.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[1], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "), 
                    # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                    "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                    # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
                    # set the migration rate between different ecotypes inhabiting the same location
                    "-m 3 1", mig_CW, "-m 4 2", mig_CW, "-m 1 3", mig_WC, "-m 2 4", mig_WC,
                    # now, set the migration rate right before the split event to zero by using the switch:
                    # -eM <t> <M>: assume a symmetric migration rate of M/(npop-1) at time t.
                    "-eM", Rsplit, "0",
                    # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                    # looking backwards in time, it moves all lines from population j into population i at time t 
                    # Migration rates into population j are set to 0 for the time further back into the past
                    "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3", 
                    # set the size of the ancient populations
                    # -en <t> <i> <n>: Set the size of population i to n*N0 at time t.
                    "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2, 
                    # add a split event - this event creates the two ancestral populations
                    "-ej", Asplit, "2 3", 
                    # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                    # -eN <t> <n>: set the size of all populations to n*N0 at time t.
                    "-eN", Asplit, 1)
  
  # create a command line for the loci without any migration (between the different ecotypes)
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and 
  # nrep: the number of independent loci that will be produced
  without.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[2], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "), "0",
                       # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                       "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                       # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                       "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3", 
                       # set the size of the ancient populations
                       # -en <t> <i> <n>: Set the size of population i to n*N0 at time t.
                       "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2, 
                       # add a split event - this event creates the two ancestral populations
                       "-ej", Asplit, "2 3", 
                       # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                       # -eN <t> <n>: set the size of all populations to n*N0 at time t.
                       "-eN", Asplit, 1)
  
  # Combine the two different types of commands
  cmdSingle <- c(with.mig, without.mig)
  cmdSingle
}

# CMD PARALLEL
cmdParallel <- function(parameters, nSites, nLoci, nDip, nPops, mutrate) {
  
  # Read the vector with the parameters and assign each parameter to the correct command name
  Ne <- parameters[1]
  # set the relative size of each population - for the extant populations
  N1 <- parameters[2]; N2 <- parameters[3]; N3 <- parameters[4]; N4 <- parameters[5]
  # and the ancient populations
  NA1 <- parameters[6]; NA2 <- parameters[7]
  
  # get the proportion of loci with migration - no barriers against migration between the different ecotypes
  pM <- parameters[14]
  # and the proportion of loci without any migration - total barrier against migration between the different ecotypes
  pNO <- parameters[15]
  
  # get the migration rates - between different ecotypes at the same site - this is the m value on the M = 4N0m formula
  # from the crab ecotype to the wave ecotype
  mCW <-  parameters[12]
  # from the wave ecotype to the crab ecotype
  mWC <-  parameters[13]
  
  # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time) 
  # and REMEMBER that M = 4N0m
  # set the migration rates to the scale of Nref - between different ecotypes at the same site
  mig_CW <- 4*Ne*mCW # -m 2 1 mig_CW
  
  # the migration from wave to crab is parametrized as a ratio of the migration from crab to wave
  # so we need to multiply the migration from crab to wave by this ratio
  mig_WC <- 4*Ne*mWC # -m 1 2 mig_WC
  
  # get the time of the recent split event
  Rsplit <- round(parameters[8], digits = 3)
  # and the ancient split event
  Asplit <- round(parameters[9], digits = 3)
  # the actual time of the ancient split event is obtained by doing:
  Asplit <- Rsplit + Asplit
  
  # Compute the value of theta
  mutrate_locus <- nSites*mutrate
  theta <- 4*Ne*mutrate_locus
  # Create a vector with information about how many haplotypes are sampled from each population
  n <- c(rep(nDip/nPops, times = nPops))
  
  # create a variable to be used for a minor correction related to the split time
  # this will ensure that the changes in migration rates occur at different times
  # correctionTime <- 1e-10
  
  # use a multinomial distribution to get the number of loci simulated under each category
  lociTotal <- as.vector(rmultinom(n = 1, size = nLoci, c(pM, pNO)))
  
  ## cheat code: pop1 - crab in site 1; pop2 - wave in site 1; pop3 - crab in site 2; pop4 - wave in site 2
  # create command line with no barriers to migration
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and 
  # nrep: the number of independent loci that will be produced
  with.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[1], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "), 
                    # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                    "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                    # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
                    # set the migration rate between different ecotypes inhabiting the same location
                    "-m 2 1", mig_CW, "-m 4 3", mig_CW, "-m 1 2", mig_WC, "-m 3 4", mig_WC,
                    # now, set the migration rate right before the split event to zero by using the switch:
                    # -eM <t> <M>: assume a symmetric migration rate of M/(npop-1) at time t.
                    "-eM", Rsplit, "0",
                    # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                    "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3", 
                    # set the size of the ancient populations
                    # -en <t> <i> <n>: Set the size of population i to n*N0 at time t.
                    "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2, 
                    # add a split event - this event creates the two ancestral populations
                    "-ej", Asplit, "2 3", 
                    # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                    # -eN <t> <n>: set the size of all populations to n*N0 at time t.
                    "-eN", Asplit, 1)
  
  # create a command line for the loci without any migration (between the different ecotypes)
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and 
  # nrep: the number of independent loci that will be produced
  without.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[2], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "), "0",
                       # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                       "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                       # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                       "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3", 
                       # set the size of the ancient populations
                       # -en <t> <i> <n>: Set the size of population i to n*N0 at time t.
                       "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2, 
                       # add a split event - this event creates the two ancestral populations
                       "-ej", Asplit, "2 3", 
                       # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                       # -eN <t> <n>: set the size of all populations to n*N0 at time t.
                       "-eN", Asplit, 1)
  
  # Combine the two different types of commands
  cmdParallel <- c(with.mig, without.mig)
  cmdParallel
}


# run SCRM --------------------------------------------------------------

# RUNSCRM
# This function will run the SCRM package and create genotypes
runSCRM <- function(Commands, nDip, nPops, model) {
  
  # check if the input is correct - the selected model should be one of the following
  if(model %in% c("Single", "Parallel", "star", "noMig", "2pops", "simple2pops") == FALSE) 
    stop(paste("The selected model should be either Single, Parallel, star, noMig, 2pops or simple2pops. Please check"))
  
  # Load the scrm package
  require(scrm)
  
  # when dealing with the model with 2 populations and no migration, the steps are more simple
  if (model == "noMig") {
    
    # Run the scrm package for each set of commands - there is only one set of commands 
    Simulation <- scrm(Commands)
    
    # Combine both simulations into a matrix of haplotypes
    Haplotypes <- unlist(Simulation, recursive = FALSE, use.names = FALSE)
    
  } else { # when dealing with the other models
    
    # Run the scrm package for each set of commands - with and without migration
    Simulation <- lapply(Commands, FUN = function(x) scrm(x))
    
    # Extract the information from each simulation and store it on a temporary matrix
    for (i in 1:length(Simulation)) {
      assign(paste("temp", i, sep = ""), Simulation[[i]][["seg_sites"]])
    }
    
    # Combine both simulations into a matrix of haplotypes
    Haplotypes <- append(temp1, unlist(mget(paste0("temp", 2:length(Simulation))), recursive = FALSE, use.names = FALSE))
  }
  
  # get the total number of haplotypes
  nHap <- nDip*2
  
  # This applies a correction for the situations where scrm does not produce a single polymorphic site
  Size <- matrix(unlist(lapply(Haplotypes, dim)), ncol=2, byrow=TRUE)
  if (length(Haplotypes[which(Size[,2]==0)]) != 0) {
    Haplotypes[which(Size[,2]==0)] <- list(matrix(rep(0, times = nHap), ncol = 1))
  }
  # It creates two sites where the haplotype is 0 for all individuals
  Size <- matrix(unlist(lapply(Haplotypes, dim)), ncol=2, byrow=TRUE)
  if (length(Haplotypes[which(Size[,2]==1)]) != 0) {
    Haplotypes[which(Size[,2]==1)] <- 
      sapply(Haplotypes[which(Size[,2]==1)], function(x) 
        list(cbind(as.matrix(unlist(x), byrow = TRUE), matrix(rep(0, times = nHap)))))
  }
  
  ## Re-organize output for the single model ##
  if (model == "Single") {
    Haplotypes <- lapply(Haplotypes, function(segSites) organizeSCRM(segSites, nHap, nPops))
  }
  
  # Convert the haplotypes to genotypes
  Genotypes <- GetGenotypes(Haplotypes, nDip = nDip)
  
  # output the genotypes
  Genotypes
}

# HAP2GENO
# Convert Haplotypes to Genotypes - this function adds the entries on one row with the entries of the subsequent row
hap2geno <- function(haplo) {
  result <- haplo[seq(1,by=2,to=nrow(haplo)),]+haplo[seq(2,by=2,to=nrow(haplo)), , drop = FALSE]
  result 
}

# GETGENOTYPES
# Create Genotypes - apply the hap2geno to all the entries of a list (different entries correspond to different loci)
GetGenotypes <- function(Haplotypes, nDip) {
  # Apply the hap2geno function across all entries of the list - sum haplotypes to get genotypes
  Genotypes <- lapply(Haplotypes, function(x) hap2geno(x))
  # Remove the name (position) of each site - this is something that scrm creates
  Genotypes <- lapply(Genotypes, function(x) unname(x))
  # Sometimes scrm creates an output with zero polymorphic sites: 
  # this adds a single column with 0s across all individuals for those situations
  Size <- matrix(unlist(lapply(Genotypes, dim)), ncol=2, byrow=TRUE)
  Genotypes[which(Size[,2]==0)] <- list(matrix(rep(0, times=nDip), ncol=1))
  Genotypes
}

# ORGANIZESCRM
# Organize SCRM output - this is utilized when dealing with the single origin model
# it switches the order of the second and third population in the final output
organizeSCRM <- function(seg_sites, nHap, nPops) {
  # Get the number of haplotypes simulated by population
  haPop <- nHap/nPops
  # Create a vector with the index representing the beginning of each population
  beginPop <- seq(from = 1, to = nHap, by = haPop)
  # Remove the name (position) of each site - this is something that scrm creates
  seg_sites <- unname(seg_sites)
  # In the single model, we need to switch the order of the second and third population
  # get the haplotypes corresponding to each population
  pop2 <- seg_sites[(beginPop[2]):(beginPop[3]-1), ]
  pop3 <- seg_sites[(beginPop[3]):(beginPop[4]-1), ]
  # re-organize the matrix of haplotypes with the populations in the correct order
  seg_sites[(beginPop[2]):(beginPop[3]-1), ] <- pop3
  seg_sites[(beginPop[3]):(beginPop[4]-1), ] <- pop2
  seg_sites
}


# calculate expected heterozygosity from genotypes -----------------------------------------------

# INDEX_INDS
# create a function to compute and store the index of individuals per population
# poolSize: a list with one entry for each population. Each entry should contain information about how many pools where used
# to sequence that population and how many individuals where in each pool.
# Example: if you used 5 pools, each with 20 individuals, the entry should be a vector such as c(20, 20, 20, 20, 20) 
index_inds <- function(poolSize) {
  
  # get the number of pops
  nPops <- length(poolSize)
  
  # set the number of individuals in the beginning - 
  # this will get updated each time the index of the individuals of a single population is computed
  nDip <- 1
  
  # create a list to store the output of the function
  sample_inds <- list()
  
  # perform a loop over all populations on the dataset
  for (i in 1:nPops) {
    
    # create the index of the individuals for any given population 
    temp <- seq.int(from = nDip, to = sum(poolSize[[i]]) + nDip - 1)
    
    # update the number of individuals - the next population should start right after the last individual of the previous one
    nDip <- nDip + max(temp)
    
    # store the index of the individuals for each population on a separated list entry
    sample_inds[[i]] <- temp
  }
  
  # output a list where each entry contains the index of the individuals for a given population
  sample_inds
}

# EXPHET_SITE
# computes the expected heterozygosity for a given site
# INPUT:
#   geno_site : vector of size nind with the genotypes coded as 0,1,2 and NA (missing data)
# OUTPUT:
#   expected heterozygosity
ExpHet_site <- function(geno_site) {  
  # get the number of individuals with data
  ngenecopies <- 2*sum(!is.na(geno_site))  
  # get the frequency of the alternative allele
  freq <- sum(geno_site, na.rm = T)/ngenecopies  
  # output the expected heterozygosity
  he <- (ngenecopies/(ngenecopies-1))*2*freq*(1-freq)
  he
}

# EXPHET
# computes the average expected heterozygosity across all sites
# INPUT:
#   geno_site : matrix of size nsites x nind with the genotypes coded as 0,1,2 and NA (missing data)
#   sample_inds: list where each entry is a vector of size nind_pop with the index of individuals belonging to a given pop
#   pop_names: vector of size npop with the names of each population
# OUTPUT:
#   average expected heterozygosity
ExpHet <- function(geno_matrix, sample_inds, pop.names) {  
  # initialize matrix to save heterozygosity for each SNP and for each pop
  # het is a matrix with nsnps rows and npop columns
  exp_het <- matrix(NA,nrow = nrow(geno_matrix), ncol = length(pop.names))
  # go through each population
  for(i in 1:length(pop.names)) {
    # compute the heterozygosity by looking at a subset of individuals from geno matrix
    exp_het[,i] <- apply(geno_matrix[,sample_inds[[i]]], 1, function(row) {ExpHet_site(row)})
  }
  exp_het
}

# EXPHET_MEAN
# computes the average expected heterozygosity across all sites
# INPUT:
#   geno_site : matrix of size nsites x nind with the genotypes coded as 0,1,2 and NA (missing data)
#   sample_inds: list where each entry is a vector of size nind_pop with the index of individuals belonging to a given pop
#   pop_names: vector of size npop with the names of each population
# OUTPUT:
#   average expected heterozygosity
ExpHet_mean <- function(geno_matrix, sample_inds, pop.names) {  
  # initialize matrix to save heterozygosity for each SNP and for each pop
  # het is a matrix with nsnps rows and npop columns
  exp_het <- matrix(NA,nrow = nrow(geno_matrix), ncol = length(pop.names))
  # go through each population
  for(i in 1:length(pop.names)) {
    # compute the heterozygosity by looking at a subset of individuals from geno matrix
    exp_het[,i] <- apply(geno_matrix[,sample_inds[[i]]], 1, function(row) {ExpHet_site(row)})
  }
  # get the mean het for each pop
  mean_exp_het <- colMeans(exp_het, na.rm = TRUE)
  names(mean_exp_het) <- pop.names
  mean_exp_het
}


# Depth of coverage -------------------------------------------------------

# SIMULATECOVERAGE
# Simulate total number of reads per site - this outputs one value per site - the depth of coverage at that site
simulateCoverage <- function(mean, variance, genotypes) {
  
  # check if the input is correct - genotypes should always be supplied as a list
  if(class(genotypes) != "list")
    stop(paste("Genotypes should be supplied on a list format, with each entry corresponding to a locus. Please check"))
  
  # Check if the variance and mean are reasonable 
  if(any(variance - mean > 0) == FALSE) {
    stop("Error: variance equal to mean, or variance smaller than mean.")
  }
  
  # calculate the parameters for the negative binomial
  pnb <- mean/variance
  rnb <- (mean^2)/(variance - mean)
  
  # Use a negative binomial to draw random values, per site and per population, for the total number of observed reads 
  # this outputs a list where each entry corresponds to a locus
  readnumbers <- lapply(genotypes, FUN = function(geno) 
    t(mapply(FUN = function(size, prob) rnbinom(n = ncol(geno), size = size, prob = prob), rnb, pnb)))
  
  # get the output - number of reads per site and per population
  readnumbers
}

# this function removes sites that have a coverage below a minimum value and sites with a coverage above a maximum value
remove_by_reads_matrix <- function(reads, minimum, maximum, genotypes = NULL) {
  
  # check which sites, if any, have a coverage below or above the threshold
  toremove <- apply(X = reads, MARGIN = 2, function(col) any(col < minimum | col > maximum))
  
  # if there are any sites with a coverage below or above the required value - remove those sites from the matrix
  if(length(toremove) != 0)
    reads <- reads[, !toremove, drop = FALSE]
  
  # set the output if the only input were the number of reads
  output <- reads
  
  # when genotypes were also supplied as input
  if(!is.null(genotypes)) {
    # remove the same sites from the matrix containing the genotypes
    genotypes <- genotypes[, !toremove, drop = FALSE]
    # create the new output, containing both the reads and the genotypes
    output <- list(reads, genotypes)
  }
  
  # output the number of reads
  output
}

# this function utilizes the previous one to remove sites above or below a threshold from multiple loci i.e. over a list
remove_by_reads <- function(nLoci, reads, minimum, maximum, genotypes = NULL) {
  
  # check if the input is correct - reads should always be supplied as a list
  if(class(reads) != "list")
    stop(paste("reads should be supplied on a list format, with each entry corresponding to a locus. Please check"))
  
  # remove sites with a depth of coverage above or below the defined threshold
  # this applies the remove_by_reads_matrix function to all the list entries - i.e. to all the different loci
  # note that this will also remove those sites from the genotypes list - if genotypes are supplied as input
  out <- lapply(1:nLoci, function(locus) 
    remove_by_reads_matrix(reads = reads[[locus]], minimum, maximum, genotypes = genotypes[[locus]]))
  
  # output the number of reads - and genotypes if relevant
  out
}

# sample total number of reads per site - this outputs one value per site - the depth of coverage at that site
sampleCoverage <- function(real, genotypes) {
  
  # check if the input is correct - real information should be on an object of class "histogram"
  if(class(real) != "histogram")
    stop("real should be an object of class histogram! 
         \nAssign the histogram of the depth of coverage of the real data to an object and use that as input")
  
  # using the information contained in the histogram:
  # sample values that are contained in the distribution of the depth of coverage from real data  
  readnumbers <- round(sample(size = ncol(genotypes), x = real$mids, prob = real$density, replace = TRUE))
  
  # get the output - number of reads per site and per population
  readnumbers 
}


# correct alpha_i ---------------------------------------------------------

set_alpha <- function(alpha_i) {
  
  # define thresholds for the minimum and maximum alpha values
  min.alpha <- 1e-2
  
  # if any value of alpha_i is below the minimum allowed value, set alpha to min.alpha value
  alpha_i[alpha_i < min.alpha] <- min.alpha
  
  # output the alpha_i for the Dirichlet distribution
  alpha_i
}


# probabilities of each pool (πp) ---------------------------------------------

# PROBSPOOL
# Calculate the probability of contribution of each pool towards the total depth of coverage of a single population
# When multiple pools were used to sequence a single population, it is possible that some pools contribute more than others
poolProbs <- function(nPools, vector_np, nSNPs, pError) {
  
  # Now, we need to calculate the probability of contributing for each individual - in each population
  # This package contains a function to perform random draws from a Dirichlet distribution
  suppressMessages(require(MCMCpack, quietly = TRUE))  
  
  # check if we are dealing with a single population - this function should be used on a single population
  # also check if the input is on the correct format
  if(class(vector_np) != "numeric" | length(vector_np) == 1)
    stop(paste("The vector_np input should be a vector. It should also contain more than one entry. Please check"))
  
  # check if we are dealing with a single population - this function should be used on a single population
  # also check if the input is on the correct format
  if(nPools != length(vector_np))
    stop(paste("The nPools input is", paste(nPools), "and so, the length of the size vector should also be", 
               paste(nPools, ".", sep = ""), "Please check"))
  
  # Silence warnings - when using only two pools, the rdirichlet functions prints a warning about reducing
  # to a beta function. This ensures that the warning is not printed
  options(warn = -1)
  
  # the total number of individuals in the population (n) can be obtained by adding the individuals in all pools
  n <- sum(vector_np)
  
  # pooling error is defined in % - change this to a proportion
  pError <- pError/100
  
  # If we use Dir(ρ*np/n), then the alpha_i for Dirichlet can be written as
  numerator <- (n - 1 - (pError^2))*vector_np
  denominator <- n*(pError^2)
  alpha <- numerator/denominator
  
  # check, and correct if needed, if any alpha_i value is above or below the threshold
  alpha <- set_alpha(alpha_i = alpha)
  
  # use a Dirichlet distribution to get the probability of contribution for each pool across all sites
  probs <- t(rdirichlet(n = nSNPs, alpha = alpha))
  
  # set the warning level back to normal
  options(warn = 0)
  
  # output the probability of contributing for each pool
  probs
}


# reads contributed by each pool (rp) -------------------------------------

# POOL READS
# Simulates the contribution, in terms of reads, for each pool 
# that contributes towards the total coverage of a given population
poolReads <- function(nPools, coverage, probs) {
  
  # set the output of this function if a given locus has no sites
  if(length(coverage) == 0) {
    
    # set the contribution to NA
    contribution <- NA
    
  } else {
    
    # simulate the contribution of each pool towards the total depth of coverage
    contribution <- vapply(1:length(coverage), FUN = function(i) {
      rmultinom(1, size = coverage[i], prob = probs[, i])
    }, FUN.VALUE = numeric(nPools))
  }
  
  # output a matrix containing the contribution (in number of reads) for each pool and across all sites
  contribution
}


# probabilities of each individual inside a pool (πip) ------------------------

# PROBSINDIVIDUAL
# Compute probability of contributing reads for each individual inside a pool
# This function works for a single pool and computes the probability of contribution for each individual inside that pool
indProbs <- function(np, nSNPs, pError) {
  
  # This package contains a function to perform random draws from a Dirichlet distribution
  suppressMessages(require(MCMCpack, quietly = TRUE))  
  
  # Silence warnings - when using only two individuals, the rdirichlet functions prints a warning about reducing
  # to a beta function. This ensures that the warning is not printed
  options(warn = -1)
  
  # pooling error is defined in % - change this to a proportion
  pError <- pError/100
  
  # If we use Dir(ρ/np), then the alpha_i for Dirichlet can be written as
  numerator <- (np - 1 - (pError^2))
  denominator <- np*(pError^2)
  alpha <- numerator/denominator
  
  # check, and correct if needed, if any alpha_i value is above or below the threshold
  alpha <- set_alpha(alpha_i = alpha)
  
  # use a dirichlet distribution to get the probability of contribution for each individual across all sites
  probs <- t(rdirichlet(n = nSNPs, alpha = rep(alpha, times = np)))
  
  # set the warning level back to normal
  options(warn = 0)
  
  # output the probability of contributing for each individual
  probs
}


# reads contributed by each individual inside a pool (rip) -------------------------------------

# INDS READS
# Simulate the contribution, in terms of reads, for each individual inside a pool
indReads <- function(nDip, coverage, probs) {
  
  # set the output of this function if a given locus has no sites
  if(length(coverage) == 0) {
    
    # set the contribution to NA
    contribution <- NA
    
  } else {
    
    # simulate the contribution of each individual towards the total depth of coverage
    contribution <- vapply(1:length(coverage), FUN = function(i) {
      rmultinom(1, size = coverage[i], prob = probs[, i])
    }, FUN.VALUE = numeric(nDip))
  }
  
  # output a matrix containing the contribution (in number of reads) for each individual and across all sites
  contribution
}


# total coverage of a given population -------------------------------------

# INDIVIDUALCOVERAGE
# compute number of reads observed for each individual and across all sites
# while the indsReads function computes the contribution of each individual inside a single pool
# this function computes the contribution of each individual inside a single population

# so, if multiple pools were used to sequence a population, this will compute the contribution of each individual inside 
# each of those pools towards the total coverage of the population

# the probability of contribution of each pool is computed and then used to calculate how many reads does that pool produce
# then the probability of contribution of each individual is computed and utilized to calculate the number of reads that 
# each individual contributes towards the total number of reads observed in the corresponding pool
popReads <- function(nPops, vector_np, coverage, pError) {
  
  # when the coverage input is a list, convert it to a vector
  if(class(coverage) == "list") 
    coverage <- unlist(coverage)
  
  # get the number of diploid individuals
  nDip <- sum(vector_np) 
  
  # get the number of pools used to sequence the population
  nPools <- length(vector_np) 
  
  # get the number of polymorphic sites 
  nSNPs <- length(coverage)
  
  # if no polymorphic sites exist at any given locus - set the output to NA
  if(nSNPs == 0) {
    # set the output 
    indCoverage <- matrix(data = NA, nrow = nPops)
    
  } else {
    
    # when the population was sequenced using a single pool of individuals
    if (nPools == 1) {
      
      # calculate the probability of contribution for each individual 
      probs <- indProbs(np = nDip, nSNPs = nSNPs, pError = pError)
      
      # compute the contribution of each individual and across all sites 
      # towards the total depth of coverage of the population
      indCoverage <- indReads(nDip = nDip, coverage = coverage, probs = probs)
      
    } else { # When more than one pool of individuals was used to sequence the population
      
      # calculate the probability of contribution for each pool used to sequence the population
      probs_pool <- poolProbs(nPools = nPools, vector_np = vector_np, nSNPs = nSNPs, pError = pError)
      
      # compute the number of reads that each pool contributes per site
      # This takes into account the total depth of coverage of the population for a given site
      # and divides that coverage amongst all the pools - considering also the probability of contribution of each pool
      pool_reads <- poolReads(nPools = nPools, coverage, probs = probs_pool)
      
      # calculate individual probabilities for each pool
      probs_ind <- lapply(vector_np, FUN = function(pool) indProbs(np = pool, nSNPs = nSNPs, pError = pError))
      
      # Taking into account the coverage of each pool, simulate how many reads does each individual inside a pool
      # contributes towards that number i.e if a pool has 25 reads at a given site: 
      # simulate how many of these reads come from the first individual, how many from the second, etc
      indCoverage <- lapply(1:nPools, function(i) 
        indReads(nDip = vector_np[i], coverage = pool_reads[i, ], probs = probs_ind[[i]]))
      
      # combine the entries from the various pools to obtain how many reads does a given individual 
      # contributes to the total depth of coverage of the population
      indCoverage <- do.call(rbind, indCoverage)
    }
  }
  
  # output a matrix containing the contribution (in number of reads) for each individual and across all sites
  indCoverage
}


# total coverage of multiple populations -------------------------------------

# POPS READS
# Compute number of reads observed for each individual and across all sites
# This function will compute the individual contribution towards the depth of coverage for each population 
# The output will be a list where each entry contains the individual contributions at each site for one of the populations
popsReads <- function(list_np, coverage, pError) {
  
  # get the number of populations
  nPops <- length(list_np)
  
  # when the coverage input is a list, convert with to a matrix, with each row corresponding to a population
  if(class(coverage) == "list") 
    coverage <-  matrix(unlist(coverage), nrow = nPops, byrow = FALSE)
  
  # Taking into account the depth of coverage for each population, simulate how many reads does each individual contributes
  # towards the total coverage of the population
  indCoverage <- lapply(1:nPops, function(pop) popReads(nPops, vector_np = list_np[[pop]], coverage[pop, ], pError))
  
  # output the contribution (in number of reads) for each individual and across all sites
  indCoverage
}


# ancestral reads ---------------------------------------------------------

# SPLITMATRIX
# Split matrix of genotypes
# This function splits a matrix into different list entries - the split is done according to row indexes 
# The size input is utilized to create the index of the rows that go into the different list entries
splitMatrix <- function(matrix, size) {
  
  # general check to see if the input is correctly supplied
  if(class(size) != "list")
    stop(paste("The size input should be a list. Please check"))
  
  # get the number of populations - it's the number of entries in the size input 
  nPops <- length(size) 
  
  # Perform a cumulative sum - this will create a vector, starting at the number one
  # each subsequent entry on the vector is the index of the last individual of a population
  popsize <- c(0, cumsum(lapply(size, sum)))
  
  # Use the previous index to create vectors containing the index of all individuals per population.
  # For example, if you have a population starting at index 30 and ending at 50, this will create a vector
  # containing the numbers 30, 31, 32... etc until 50. This creates a list, where each entry is a vector with  
  # the indices of a single population
  index <- lapply(1:nPops, function(i) seq(popsize[i]+1, popsize[i+1]))
  
  # use the vectors containing the index of all individuals of a given population to subset the matrix
  # The subsetting is done by rows and it creates a list, where each entry contains the information for one population
  output <- lapply(1:nPops, function(pop) matrix[index[[pop]], , drop = FALSE])
  
  # Output the list containing, in each entry, the information for each population
  output
}  

# GETNUMREADSA_VECTOR
# Calculate the number of Ancestral reads
# this function computes the number of reads with the ancestral allele over a vector containing depths of coverage
getNumReadsA_vector <- function(genotype_v, readCount_v, error) {
  
  # ensure that both vectors have the same length
  # you should have one genotype and one value of depth of coverage per site
  if(length(genotype_v) != length(readCount_v))
    stop(paste("The lengths of the genotype_v and the readCount_v inputs should be the same. Please check"))
  
  # initialize the number of A reads as a vector of 0 with same size as Genotype_v
  numAreads <- numeric(length(genotype_v))
  
  # Get the entries of genotype for hom ancestra
  hom_anc <- genotype_v == 0
  # Get the entries of heterozygotes (1)
  het <- genotype_v == 1
  # Get the entries for hom derived (2)
  hom_der <- genotype_v == 2
  
  # Call the rbinom in a vector way, using size as the vector
  # homozygote ancestral
  numAreads[hom_anc] <- rbinom(n = sum(hom_anc), size = readCount_v[hom_anc], prob = 1-error)
  # hets
  numAreads[het] <- rbinom(n = sum(het), size = readCount_v[het], prob = 0.5)
  # homozygote derived
  numAreads[hom_der] <- rbinom(n = sum(hom_der), size = readCount_v[hom_der], prob = error)
  # return the number of reads
  numAreads
}

# COMPUTEANCESTRAL
# Calculate the number of Ancestral reads over a matrix
# Both the genotypes and the indContribution inputs are matrices - and they should have the same dimensions
# this function treats each row of the matrices as a vector and computes the number of reads with the ancestral allele 
computeAncestral <- function(genotypes, indContribution, error) {
  
  # # ensure that both matrices have the same dimension
  # you should have one genotype and one value of depth of coverage per site - both matrices should have the same dimension
  if(identical(dim(genotypes), dim(indContribution)) == FALSE)
    stop(paste("The dimensions of the genotypes and indContribution matrices should be the same. Please check"))
  
  # Vectorize the function getNumReadsA_vector
  # The function getNumReadsA_vector gets as input vectors
  # Hence, I can apply it to each column of genotypes and individual contribution
  tempAncestral <- sapply(1:ncol(genotypes), function(i) {
    getNumReadsA_vector(genotype_v = genotypes[,i], readCount_v = indContribution[,i], error = error)})
  
  # Check that indContribution-tmp does not give negative values
  if(sum( (indContribution - tempAncestral) < 0 ) > 0 ) {
    stop("Error in creating number of reads!")
  }
  
  # output the number of ancestral reads
  tempAncestral
}

# NUMBERANCESTRAL
# Calculate the number of Ancestral reads over a list - where each entry of the genotypes list is a different locus
# this also works when using a single locus - provided that the input is in the list format
# note that each entry of the genotypes list should be a matrix of genotypes  
# IMPORTANT: this function is coded to work on a single population!
numberAncestral <- function(genotypes, indContribution, error) {
  
  if(class(genotypes) != "list")
    stop(paste("The genotypes input should be on a list format. Please check"))
  
  if(class(indContribution) != "list")
    indContribution <- list(indContribution)
  
  if(length(genotypes) != length(indContribution))
    stop(paste("The genotypes and indContribution lists should have the same number of entries. Please check"))
  
  # get the number of loci 
  nLoci <- length(genotypes)
  
  # When dealing with a single population - all the genotypes in the matrix correspond to a single population
  readsAncestral <- lapply(1:nLoci, function(locus) computeAncestral(genotypes[[locus]], indContribution[[locus]], error))
  
  # output the number of reads with the ancestral allele per individual and per site
  readsAncestral
} 

# NUMBERANCESTRALPOP
# Calculate the number of Ancestral reads for multiple populations
# genotypes could be either a matrix of genotypes - coded as 0s, 1s and 2s 
# and where each column is a site and each row one individual
# OR a list with a single entry (one locus) - that entry should be a matrix of genotypes
# indContribution is a list where each entry contains information for a single population
# each entry of that list should be a matrix with the number of reads observed for each individual
# IMPORTANT: this function is coded to work on more than one population but for a single locus!
numberAncestralPop <- function(genotypes, indContribution, size, error) {
  
  # get the number of populations
  nPops <- length(size)
  
  # check if the indContribution input is the right format
  if(class(indContribution) != "list")
    stop(paste("The indContribution input should be on a list format. Please check"))
  
  # set the output when a given locus has no SNPs
  if(any(is.na(unlist(indContribution)))) {
    # set the output to NA
    ancestral <- list(matrix(data = NA, nrow = nPops), matrix(data = NA, nrow = nPops))
    
  } else {
    
    # when the genotypes input is not a list, convert with to a list
    if(class(genotypes) != "list") 
      genotypes <-  list(genotypes)
    
    # check if the indContribution list contains one list entry per population
    if(length(indContribution) != nPops)
      stop(paste("The indContribution input should have one list entry per population. Please check"))
    
    # use the splitMatrix function to split the genotypes into different list entries for each population
    tempGeno <- sapply(genotypes, FUN = function(geno) splitMatrix(geno, size))
    
    # compute the number of reads with the ancestral allele, for each individual and across all sites and populations
    ancestral <- sapply(1:nPops, function(pop)
      numberAncestral(genotypes = list(tempGeno[[pop]]), indContribution = indContribution[[pop]], error))
  }
  
  # output the number of ancestral reads
  ancestral
}


# create pools ------------------------------------------------------------

# POOL_ONEPOP
# Create Pooled DNA sequencing data
# this function combines the information for each individual - contained in the indContribution and readsAncestral inputs
# into information at the population level - i.e. the information of all individuals in the population is combined 
# into a single population value 
pool_onePop <- function(nLoci, indContribution, readsAncestral) {
  
  # when dealing with a single population and one locus it's possible that the indContribution input is not on a list format
  # if this is the case, transform it into a list
  if(class(indContribution) != "list")
    indContribution <- list(indContribution) 
  
  # Number of reads with the derived allele is simply the total number of reads per individual minus the number of reads 
  # with the ancestral allele for that individual
  readsDerived <- lapply(1:nLoci, function(locus) indContribution[[locus]] - readsAncestral[[locus]])
  
  # Create matrices to store the number of total reads and reads with ancestral and derived alleles per pool
  # note that nrow = 1 - because we are only working with one population
  ancestralPool <- lapply(readsAncestral, function(x) matrix(colSums(x), nrow = 1))
  derivedPool <- lapply(readsDerived, function(x) matrix(colSums(x), nrow = 1))
  totalPool <- lapply(indContribution, function(x) matrix(colSums(x), nrow = 1))
  
  # Combine the information about the different read types into a list, create names for each entry and output that list
  final_list <- list(ancestralPool, derivedPool, totalPool)
  names(final_list)=c("Ancestral", "Derived", "Total")
  return(final_list)
}

# POOLPOPS
# The following is used when working with more than one population
# in this situation, each entry of the indContribution and readsAncestral lists should contain one entry per population
# so, it is, in essence, a list within a list
# this function then combines the individual information into population information
poolPops <- function(nPops, nLoci, indContribution, readsAncestral) {
  
  # when dealing with a single locus it's possible that the indContribution input is not on a list format
  if(nLoci == 1) {
    indContribution <- list(indContribution); readsAncestral <- list(readsAncestral)
  }
  
  # Number of reads with the derived allele is simply the total number of reads per individual minus the number of reads 
  # with the ancestral allele for that individual
  readsDerived <- lapply(1:nLoci, function(locus) 
    mapply(function(individual, ancestral) FUN = individual - ancestral, SIMPLIFY = FALSE, 
           individual = indContribution[[locus]], ancestral = readsAncestral[[locus]]))
  
  # Now, since each entry (each locus) has independent entries for each population, we can simple perform colSums across the 
  # various entries. This will sum the number of reads (with the ancestral, the derived and the total number of reads) 
  # across all individuals of a population.
  # in this way, you get the total number of reads with each allele for each population
  ancestralPool <- lapply(readsAncestral, function(ancestral) matrix(t(sapply(ancestral, colSums)), nrow = nPops))
  derivedPool <- lapply(readsDerived, function(derived) matrix(t(sapply(derived, colSums)), nrow = nPops))
  totalPool <- lapply(indContribution, function(total) matrix(t(sapply(total, colSums)), nrow = nPops))
  
  # Combine the information about the different read types into a list, create names for each entry and output that list
  final_list <- list(ancestralPool, derivedPool, totalPool)
  names(final_list)=c("Ancestral", "Derived", "Total")
  return(final_list)
}


# filter SNPs -------------------------------------------------------------

# MAJOR SIMS 
# this function is utilized to ensure that the ancestral allele of the simulations is also the major allele 
# whenever the number of derived allele reads - across all populations - is larger than the number of ancestral allele reads
# it swaps the columns of the two matrices (ancestral and derived), so that the most frequent allele is in the ancestral matrix
# then, if min.minor is not NA: 
# this function also removes sites where the number of minor allele reads - across all populations - is below a certain threshold
majorSims <- function(ancestral, derived, coverage, min.minor = NA) {
  
  # set the output for the situations where there is no SNP at the locus
  # check for NAs in one of the matrices
  if(any(is.na(unlist(ancestral)))) {
    # set the output to NA
    out <- list(Ancestral = NA, Derived = NA, Total = NA)
    
  } else {
    
    # create two temporary matrices - one containing the reads with the ancestral allele for each population and across all sites
    # and another with the reads with the derived allele
    Ranc <- ancestral; Rder <- derived
    
    # perform an evaluation to check if the total number or reads for each site 
    # is bigger in the matrix containing the "ancestral" allele reads or in the matrix containing the "derived" allele
    eval <- colSums(Ranc) < colSums(Rder)
    
    # at the sites where the number of reads for the "derived" allele is bigger than the number of reads for the "ancestral" allele
    # replace those columns with the rows from the matrix of the "derived" allele
    ancestral[, eval] <- Rder[, eval]
    # and then replace those same columns in the matrix with the reads of the "derived" allele
    derived[, eval] <- Ranc[, eval]
    
    # if the min.minor input is not NA - then it should be an integer 
    # representing the minimum number of reads with the minor allele (the derived one) that we should observe across all populations
    if (is.na(min.minor) == FALSE) {
      # now we need to find the total coverage - across all populations - of the minor frequency allele
      # since we switched the columns in the previous step, the number of reads with the minor allele - the less frequent allele
      # are stored in the derived matrix
      minor <- colSums(derived)
      # find out in which columns the total sum of the reads with the minor allele is below the threshold
      toremove <- minor < min.minor
      
      # if there are sites where the sum of the reads with the minor allele is below the threshold
      if (length(toremove) != 0) {
        # remove those columns from the matrix containing the depth of coverage
        coverage <- coverage[, !toremove, drop = FALSE]
        # remove those columns from the matrix containing the number of reads with the ancestral allele
        ancestral <- ancestral[, !toremove, drop = FALSE]
        # remove those columns from the matrix containing the number of reads with the derived allele
        derived <- derived[, !toremove, drop = FALSE]
      }
    }
    
    # output the matrices with the various types of read numbers
    out <- list(Ancestral = ancestral, Derived = derived, Total = coverage)
  }
  
  # output the matrices with the various types of read numbers
  out
}


# pick loci ---------------------------------------------------------------

# change the pick loci function to randomly select loci from the lists of the pooled data
# keeping the proportion of loci simulated with or without migration
pickLoci <- function(commands, nPops, nLoci, extra, major, minor, coverage) {
  
  # get the number of loci simulated with normal migration
  nMig <- as.numeric(sapply(strsplit(commands[1], " "), "[[", 2))
  # and the number simulated without migration
  nPNo <- as.numeric(sapply(strsplit(commands[2], " "), "[[", 2))
  # check if they sum to the expected value
  if((nMig + nPNo) != extra)
    stop("There is a mistake in the number of simulated loci")
  
  # set the object containing the coverage
  cov <- coverage
  
  # set the default value of the expected number of loci without migration:
  exp.pNo <- 0  
  # and create lists for the number of reads with the major allele, the number of reads with the minor allele 
  # and the total coverage of the loci without migration
  pNo.maj <- list(); pNo.min <- list(); pNo.cov <- list()
  
  # if the number of loci without migration is not zero
  if(nPNo != 0) {
    
    # get the proportion of the loci with or without migration
    probs <- c(nMig/extra, nPNo/extra)
    # given this proportion, in "nLoci" we would expect the number of loci without migration to be:
    exp.pNo <- round(probs[2] * nLoci)
    
    # get the major allele reads corresponding to the loci without migration
    pNo.maj <- tail(major, n = nPNo)
    # get the minor allele reads corresponding to the loci without migration
    pNo.min <- tail(minor, n = nPNo)
    # and the coverage of the loci without migration - from the coverage list
    pNo.cov <- tail(cov, n = nPNo)
    
    # get the major allele reads corresponding to the loci with migration
    major <- head(major, -nPNo)
    # get the minor allele reads corresponding to the loci with migration
    minor <- head(minor, -nPNo)
    # and the coverage of the loci with migration
    cov <- head(cov, -nPNo)
    
    # remove any empty list elements - from both the reads list and the genotypes list
    # for both sets of lists:
    # with migration between the ecotypes
    major <- major[!sapply(major, length) == 0]; minor <- minor[!sapply(minor, length) == 0] 
    cov <- cov[!sapply(cov, length) == 0]
    # and without migration between the ecotypes
    pNo.maj <- pNo.maj[!sapply(pNo.maj, length) == 0]; pNo.min <- pNo.min[!sapply(pNo.min, length) == 0]
    pNo.cov <- pNo.cov[!sapply(pNo.cov, length) == 0]
    
    # we should also remove any elements that have an NA - corresponding to locus without a SNP prior to the filtering
    # from the lists containing the information about the loci with migration
    major <- major[!is.na(major)]; minor <- minor[!is.na(minor)]; cov <- cov[!is.na(cov)]
    # and without migration between the ecotypes
    pNo.maj <- pNo.maj[!is.na(pNo.maj)]; pNo.min <- pNo.min[!is.na(pNo.min)]; pNo.cov <- pNo.cov[!is.na(pNo.cov)]
    
    # ensure that each entry is a matrix - for the loci with migration
    major <- lapply(major, function(locus) matrix(locus, nrow = nPops))
    # and without migration
    pNo.maj <- lapply(pNo.maj, function(locus) matrix(locus, nrow = nPops))
    
    # if there are not enough loci left to pick
    if(length(pNo.maj) < exp.pNo) {
      
      # set the output to NA
      pNo.maj <- NA; pNo.min <- NA; pNo.cov <- NA
      
    } else {
      
      # now, from the list containing the reads for the loci without migration
      # select random entries to keep - note that we are only keeping the number of entries expected for nLoci locus
      tokeep <- sort(sample(x = 1:length(pNo.maj), size = exp.pNo))
      # then, keep only those entries in the list containing the major allele reads without migration
      pNo.maj <- pNo.maj[tokeep]
      # and from the list containing the minor allele reads of loci without migration 
      pNo.min <- pNo.min[tokeep]
      # and in the list containing the total coverage of the loci without migration
      pNo.cov <- pNo.cov[tokeep]  
    }
  }
  
  # remove any empty list elements - from both the reads list and the coverages list
  major <- major[!sapply(major, length) == 0]; minor <- minor[!sapply(minor, length) == 0]; cov <- cov[!sapply(cov, length) == 0]
  # and any list element containing an NA
  major <- major[!is.na(major)]; minor <- minor[!is.na(minor)]; cov <- cov[!is.na(cov)]
  
  # if there are not enough loci left to pick
  if(length(major) < (nLoci - exp.pNo)) {
    
    # set the output to NA
    major <- NA; minor <- NA; cov <- NA
    
  } else {
    
    # then we need to select random loci with migration between the ecotypes to keep 
    # select random loci to keep in order to obtain a sample of nLoci size
    tokeep <- sort(sample(x = 1:length(major), size = nLoci - exp.pNo))
    
    # keep only those loci in the major allele reads list
    major <- major[tokeep]
    # in the minor allele list
    minor <- minor[tokeep]
    # and in the list with the total depth of coverage
    cov <- cov[tokeep]
  }
  
  # now, combine both the major allele reads of loci with or without migration into a single list
  major <- c(major, pNo.maj)
  # do the same for the minor allele reads
  minor <- c(minor, pNo.min)
  # and for the total coverage of loci with or without migration
  cov <- c(cov, pNo.cov)
  
  # output both the reads and the genotypes
  out <- list(ancestral = major, derived = minor, total = cov)
}

  		  
# fraction of sites -------------------------------------------------------

# this function will compute the fraction of sites showing a fixed difference between populations
  # ancestral is a matrix with the number of reads for the ancestral allele
  # total is a matrix with the total number of reads
fixed <- function(ancestral, total, nPops) {
  
  # check if the input is in the correct format
  if(nrow(ancestral) != nPops | nrow(total) != nPops)
    stop(paste("Using an incorrect input. There should be one row per population in both the ancestral and total matrices"))
  
  # get the number of fixed sites in the first population that, at the same time, do not exist in the second population
  fixA <- ancestral[1, ] == total[1, ] & ancestral[2, ] == 0
  # get the number of fixed sites in the second population that, at the same time, do not exist in the first population
  fixB <- ancestral[2, ] == total[2, ] & ancestral[1, ] == 0  
  
  # combine the previous information from both populations
  fixAB <- c(fixA, fixB)
  
  if (nPops == 2) {
    
    # fraction of sites showing a fixed difference between both populations is obtained by performing a sum of the previous vector
    # because each time one populations has a frequency of 1 and the other has a frequency of 0 - that corresponds to a TRUE
    # and dividing that sum by the total number of sites
    Sf <- sum(fixAB) / ncol(ancestral)
  
  } else {
    
    # get the number of fixed sites in the third population that, at the same time, do not exist in the fourth population
    fixC <- ancestral[3, ] == total[3, ] & ancestral[4, ] == 0
    # get the number of fixed sites in the fourth population that, at the same time, do not exist in the third population
    fixD <- ancestral[4, ] == total[4, ] & ancestral[3, ] == 0
    # combine the previous information from both populations - this is the information for the other location
    fixCD <- c(fixC, fixD)
    # compute the fraction of sites with a fixed difference between the populations at each location 
    Sf <- c(sum(fixAB) / ncol(ancestral), sum(fixCD) / ncol(ancestral))
    
    # We can also look at this globally - for all the four populations 
    # check, for each population, which sites are fixed for that population and absent from the others
    fixA <- ancestral[1, ] == total[1, ] & ancestral[2, ] == 0 & ancestral[3, ] == 0 & ancestral[4, ] == 0 
    fixB <- ancestral[2, ] == total[2, ] & ancestral[1, ] == 0 & ancestral[3, ] == 0 & ancestral[4, ] == 0 
    fixC <- ancestral[3, ] == total[3, ] & ancestral[1, ] == 0 & ancestral[2, ] == 0 & ancestral[4, ] == 0 
    fixD <- ancestral[4, ] == total[4, ] & ancestral[1, ] == 0 & ancestral[2, ] == 0 & ancestral[3, ] == 0
    # combine all of this information into a single vector
    gfix <- c(fixA, fixB, fixC, fixD)
    # the global fraction of sites with a fixed difference can be obtained by dividing the sum of the previous vector
    # by the total number of sites
    gSf <- sum(gfix) / ncol(ancestral)
    # add this global fraction to the vector containing the information about the pairwise comparison of locations
    Sf <- c(Sf, gSf)
  }
  
  # output the fraction of sites showing a fixed difference between the populations
  Sf
}

# this function will compute the fraction of sites showing an exclusive polymorphism to a given population
  # ancestral is a matrix with the number of reads for the ancestral allele
  # total is a matrix with the total number of reads
exclusive <- function(ancestral, total, nPops) {
  
  # check if the input is in the correct format
  if(nrow(ancestral) != nPops | nrow(total) != nPops)
    stop(paste("Using an incorrect input. There should be one row per population in both the ancestral and total matrices"))
  
  # get the number of segregating sites in the first population
  # i.e. the number of sites where the number of reads with the ancestral allele is not equal to either the total nreads or zero 
  segA <- ancestral[1, ] != total[1, ] & ancestral[1, ] != 0
  # get the number of segregating sites in the second population
  segB <- ancestral[2, ] != total[2, ] & ancestral[2, ] != 0
  
  # combine the previous information from both populations - this is the information for one of the locations
  segAB <- rbind(segA, segB)
  
  # we are only interested in the sites where the sum of values per column is equal to 1
  # because this means that we have TRUE (i.e segregating site) in one population and FALSE in the other (non-segregating site)
  keep <- colSums(segAB) == 1
  # remove columns where both populations are fixed for one allele or have no reads for that allele
  segAB <- segAB[, keep, drop = FALSE]
  # now by performing a sum over the rows we can obtain the number of sites, for each of the populations 
  # where that population is polymorphic while the other is not
  segAB <- unname(rowSums(segAB))

  if (nPops == 2) {
    
    # now we can calculate the fraction of sites showing an exclusive polymorphism for each population
    # this is obtained by dividing the number of exclusive sites by the total number of sites
    Sx <- segAB / ncol(ancestral)
    
  } else {
  
    # now we can calculate the fraction of sites showing an exclusive polymorphism for each population
    # this is obtained by dividing the number of exclusive sites by the total number of sites
    SxAB <- segAB / ncol(ancestral)
    
    # get the number of segregating sites in the third population
    segC <- ancestral[3, ] != total[3, ] & ancestral[3, ] != 0
    # get the number of segregating sites in the fourth population
    segD <- ancestral[4, ] != total[4, ] & ancestral[4, ] != 0
    # combine the previous information from both populations - this is the information for the second location
    segCD <- rbind(segC, segD)
    
    # we are only interested in the sites where the sum of values per column is equal to 1
    # because this means that we have TRUE (i.e segregating site) in one population and FALSE in the other (non-segregating site)
    keep <- colSums(segCD) == 1
    # remove columns where both populations are fixed for one allele or have no reads for that allele
    segCD <- segCD[, keep, drop = FALSE]
    # now by performing a sum over the rows we can obtain the number of sites, for each of the populations 
    # where that population is polymorphic while the other is not
    segCD <- unname(rowSums(segCD))
    # now we can calculate the fraction of sites showing an exclusive polymorphism for each population
    # this is obtained by dividing the number of exclusive sites by the total number of sites
    SxCD <- segCD / ncol(ancestral)
    
    # We can also look at this globally - for all the four populations 
    seg <- rbind(segA, segB, segC, segD)
    # again, we are only interested in sites where the sum of the columns is one
    # this time, this means that the allele is segregating in only one population and not in the other three
    keep <- colSums(seg) == 1
    # now by performing a sum we can obtain the number of sites where one population is polymorphic while the others are not
    seg <- sum(keep)
    # calculate the fraction of sites showing an exclusive polymorphism for a single population
    gSeg <- seg / ncol(ancestral)
    
    # combine all the different information into a single vector
    Sx <- c(SxAB, SxCD, gSeg)
  }
  
  # output the fraction of sites showing an exclusive polymorphism for a given population
  Sx
}

# this function will compute the fraction of sites with a shared polymorphism between the populations
  # ancestral is a matrix with the number of reads for the ancestral allele
  # total is a matrix with the total number of reads
shared <- function(ancestral, total, nPops) {
  
  # check if the input is in the correct format
  if(nrow(ancestral) != nPops | nrow(total) != nPops)
    stop(paste("Using an incorrect input. There should be one row per population in both the ancestral and total matrices"))
  
  # get the number of segregating sites in the first population
  # i.e. the number of sites where the number of reads with the ancestral allele is not equal to either the total nreads or zero 
  segA <- ancestral[1, ] != total[1, ] & ancestral[1, ] != 0
  # get the number of segregating sites in the second population
  segB <- ancestral[2, ] != total[2, ] & ancestral[2, ] != 0
    
  # combine the information into a single matrix
  segAB <- rbind(segA, segB)
  # we are only interested in the sites where the sum of values per column is equal to 2
  # because this means that we have TRUE (i.e segregating site) in one population and TRUE in the other (also segregating)
  sharedAB <- sum(colSums(segAB) == 2)
  # the total number of TRUES in the previous line is the number of sites with a shared polymorphism between the populations
  
  if (nPops == 2) {
    
    # calculate the fraction of sites showing a shared polymorphism between the two populations
    # this is obtained by dividing the number of shared sites by the total number of sites
    SS <- sharedAB / ncol(ancestral)
    
  } else {
    
    # get the number of segregating sites in the third population
    segC <- ancestral[3, ] != total[3, ] & ancestral[3, ] != 0
    # get the number of segregating sites in the fourth population
    segD <- ancestral[4, ] != total[4, ] & ancestral[4, ] != 0
    # combine the previous information from both populations - this is the information for the other location
    segCD <- rbind(segC, segD)
    # get number of sites with a shared polymorphism between the populations in the second location
    sharedCD <- sum(colSums(segCD) == 2)
    
    # now we can calculate the fraction of sites showing a shared polymorphism between the populations - and for each site
    # this is obtained by dividing the number of exclusive sites by the total number of sites
    SS <- c(sharedAB / ncol(ancestral), sharedCD / ncol(ancestral))
    
    # We can also look at this globally - for all the four populations 
    seg <- rbind(segA, segB, segC, segD)
    # get number of sites with a shared polymorphism between all of the populations
    shared <- sum(colSums(seg) == nPops)
    # the fraction of sites with a polymorphism shared between all of the populations is obtained by dividing
    # the number of shared sites by the total number of sites
    gSeg <- shared / ncol(ancestral)
    
    # combine all the different information into a single vector
    SS <- c(SS, gSeg)
  }   
  
  # output the fraction of sites with a polymorphism shared between the populations
  SS
}


# allelic frequencies -----------------------------------------------------

# CALCULATEPI
# Calculate population frequency at a SNP locus
# this function computes the allelic frequency - this is calculated as (nreads ancestral)/(nreads total)
# it also removes sites without information - sites where there is an NA?
calculatePi <- function(listPool, nLoci) {
  
  # The input should be the list containing information about reads per pool that was created by the previous function
  # A failsafe to detect if you're using the right list
  if (("ancestral" %in% names(listPool) | "total" %in% names(listPool)) == FALSE) {
    stop(paste("Using an incorrect list"))
  }
  
  # It is possible that each entry of the list is a single matrix - particularly when dealing with a single locus
  # If this is the case, then this step will convert those entries into lists
  if(any(lapply(listPool, class) == "matrix") == TRUE)
    listPool <- lapply(listPool, list)
  
  # By doing that transformation, we can use an lapply, whether we have one locus or multiple loci
  # Divide the number of reads with the ancestral allele by the total number of reads - to obtain allelic frequencies
  Site_Pi <- lapply(1:nLoci, function(locus) listPool[["ancestral"]][[locus]]/listPool[["total"]][[locus]])
  
  # Remove sites without information from the various categories of reads
  listPool[["ancestral"]] <- lapply(1:nLoci, function(locus) 
    listPool[["ancestral"]][[locus]][, colSums(is.na(Site_Pi[[locus]])) == 0, drop = FALSE])
  
  listPool[["derived"]] <- lapply(1:nLoci, function(locus) 
    listPool[["derived"]][[locus]][, colSums(is.na(Site_Pi[[locus]])) == 0, drop = FALSE])
  
  listPool[["total"]] <- lapply(1:nLoci, function(locus) 
    listPool[["total"]][[locus]][, colSums(is.na(Site_Pi[[locus]])) == 0, drop = FALSE])
  
  # Remove sites without information from the matrix containing the allelic frequencies
  Site_Pi <- lapply(1:nLoci, function(locus) 
    Site_Pi[[locus]][, colSums(is.na(Site_Pi[[locus]])) == 0, drop = FALSE])
  
  # The output is a list that contains the information about the site frequencies
  # But also information about the total number of reads and the number of reads with the ancestral/derived alleles 
  list(Site_Pi, listPool)
}


# mean expected heterozygosity --------------------------------------------

# EXPHET_SITE
# computes the expected heterozygosity for a given site
# INPUT:
#   geno_site : vector of size nind with the genotypes coded as 0,1,2 and NA (missing data)
# OUTPUT:
#   expected heterozygosity
ExpHet_site <- function(geno_site) {  
  # get the number of individuals with data
  ngenecopies <- 2*sum(!is.na(geno_site))  
  # get the frequency of the alternative allele
  freq <- sum(geno_site, na.rm = T)/ngenecopies  
  # output the expected heterozygosity
  he <- (ngenecopies/(ngenecopies-1))*2*freq*(1-freq)
  he
}

# EXPHET
# computes the average expected heterozygosity across all sites
# INPUT:
#   geno_site : matrix of size nsites x nind with the genotypes coded as 0,1,2 and NA (missing data)
#   sample_inds: list where each entry is a vector of size nind_pop with the index of individuals belonging to a given pop
#   pop_names: vector of size npop with the names of each population
# OUTPUT:
#   average expected heterozygosity
ExpHet <- function(geno_matrix, sample_inds, pop.names) {  
  # initialize matrix to save heterozygosity for each SNP and for each pop
  # het is a matrix with nsnps rows and npop columns
  exp_het <- matrix(NA,nrow = nrow(geno_matrix), ncol = length(pop.names))
  # go through each population
  for(i in 1:length(pop.names)) {
    # compute the heterozygosity by looking at a subset of individuals from geno matrix
    exp_het[,i] <- apply(geno_matrix[,sample_inds[[i]]], 1, function(row) {ExpHet_site(row)})
  }
  # check the output
  #str(exp_het)  
  exp_het
}

# EXPHET_MEAN
# computes the average expected heterozygosity across all sites
# INPUT:
#   geno_site : matrix of size nsites x nind with the genotypes coded as 0,1,2 and NA (missing data)
#   sample_inds: list where each entry is a vector of size nind_pop with the index of individuals belonging to a given pop
#   pop_names: vector of size npop with the names of each population
# OUTPUT:
#   average expected heterozygosity
ExpHet_mean <- function(geno_matrix, sample_inds, pop.names) {  
  # initialize matrix to save heterozygosity for each SNP and for each pop
  # het is a matrix with nsnps rows and npop columns
  exp_het <- matrix(NA,nrow = nrow(geno_matrix), ncol = length(pop.names))
  # go through each population
  for(i in 1:length(pop.names)) {
    # compute the heterozygosity by looking at a subset of individuals from geno matrix
    exp_het[,i] <- apply(geno_matrix[,sample_inds[[i]]], 1, function(row) {ExpHet_site(row)})
  }
  # check the output
  #str(exp_het)  
  # get the mean het for each pop
  mean_exp_het <- colMeans(exp_het, na.rm = TRUE)
  names(mean_exp_het) <- pop.names
  mean_exp_het
}

# NOTE: the previous functions are utilized to compute expected heterozygosity directly from genotypes
# the following functions compute expected heterozygosity from pooled data

# EXPECTED_HET
# expected heterozygosity within a population 
Expected_Het <- function(Pop_Pi) {
  # Dealing with a single matrix of population allelic frequencies - a single locus or simulation
  if(class(Pop_Pi) == "matrix") {
    # Compute the expected heterozygosity for a site - this code goes across all sites
    het <- apply(Pop_Pi, c(1,2), function (frequency) 2*frequency*(1 - frequency))
    
  } else { # Dealing with more than one locus or simulation
    
    het <- lapply (Pop_Pi, FUN = function(x) {
      apply(x, c(1,2), function (frequency) 2*frequency*(1 - frequency))})
  }
  # output the expected heterozygosity 
  het
}

# MEANEXPECTED_HET
meanExpected_Het <- function(Pop_Pi) {
  # Dealing with a single matrix of population allelic frequencies - a single locus or simulation
  if(class(Pop_Pi) == "matrix") {
    # Compute the expected heterozygosity for a site - this code goes across all sites
    het <- apply(Pop_Pi, c(1,2), function (frequency) 2*frequency*(1 - frequency))
    # Compute the mean across rows - each pool (or population) occupies one row 
    # so this gives the mean expected heterozygosity of each population
    mean_het <- rowMeans(het)
    
  } else { # Dealing with more than one locus or simulation
    
    het <- lapply (Pop_Pi, FUN = function(x) {
      apply(x, c(1,2), function (frequency) 2*frequency*(1 - frequency))})
    mean_het <- lapply (het, FUN = function(x) {rowMeans(x)})
  }
  # output the mean expected heterozygosity 
  mean_het
}


# expected heterozygosity between population pairs --------------------------------------------

# HET_BETWEEN
# this function computes the heterozygosity between populations 
# the output is a single value for each pairwise comparison
Het_Between <- function(Pop_Pi) {
  ## when dealing with a single matrix of population frequencies - one locus or one simulation
  ## #########################################################################################
  if(class(Pop_Pi) == "matrix") {
    # Considering only two populations (according to Hudson (1992) estimator, using the formula of Chen (2015):
    # http://journals.plos.org/plosone/article/metrics?id=10.1371/journal.pone.0135368
    if(nrow(Pop_Pi) == 2) {
      HB <- (1/2)*((2*Pop_Pi[1,]*(1-Pop_Pi[2,]))+(2*Pop_Pi[2,]*(1-Pop_Pi[1,])))
      # (1/2) because we are comparing pairs of populations
      Mean_HB <- mean(HB)
      
    } else { # Considering more than two populations
      Pairwise <- combn(nrow(Pop_Pi),2)
      HB <- as.matrix(apply(Pairwise, 2, function(col) {
        (1/2)*(2*Pop_Pi[col[1],]*(1-Pop_Pi[col[2],])+2*Pop_Pi[col[2],]*(1-Pop_Pi[col[1],]))}))
      
      # This adds a correction for the situations where there is only one site 
      if(ncol(HB) == 1) {
        HB <- t(HB)
      }
      
      # Compute the mean per column - each column contains the values for one of the pairwise combinations
      # So by computing the column mean, we compute the mean heterozygosity between two populations
      Mean_HB <- colMeans(HB)
    }
    
    ## Starting here, we are dealing with the situations where we have more than one simulation or locus
    ## #################################################################################################
    # Meaning that the input is a list, instead of a matrix
  } else {
    # get the number of pops - look at the number of the rows for the first simulation or locus
    nPops <- nrow(Pop_Pi[[1]])
    # create a matrix with all the possible pairwise comparisons between the populations
    Pairwise <- combn(nPops, 2)
    HB <- lapply(Pop_Pi, FUN = function(y) {
      as.matrix(apply(Pairwise, 2, function(col) {
        (1/2)*(2*y[col[1],]*(1-y[col[2],])+2*y[col[2],]*(1-y[col[1],]))}))
    })
    
    # when dealing with just one comparison - two populations - you only need to compute the mean across all sites
    if (nPops == 2) {
      Mean_HB <- lapply(HB, FUN = function(x) {mean(x)}) 
      
    } else { # when dealing with more than two populations
      
      # Correction for simulations with a single polymorphic site
      HB <- lapply(HB, FUN = function(H) {
        if(ncol(H) == 1) {
          H <- t(H)
        } else {
          H = H
        }
      })
      # Compute the mean for each column 
      Mean_HB <- lapply(HB, FUN = function(x) {colMeans(x)}) 
    }
  } 
  # output the results of the function
  Mean_HB
}


# compute FST values ------------------------------------------------------

# GETFST
# Calculate FST Values
# This function outputs a single value - FST between two populations at a single site
getFst <- function(freq1, freq2, ss1, ss2) {
  p1 <- freq1*(1-freq1)
  p2 <- freq2*(1-freq2)
  pdiffsquare <- (freq1-freq2)^2
  numerator <- pdiffsquare - (p1/(ss1-1)) - (p2/(ss2-1))
  denominator <- (freq1*(1-freq2)) + (freq2*(1-freq1))
  den <- sum(denominator, na.rm = T)
  if (den > 0) {
    res <- sum(numerator, na.rm = T)/den
  } else {
    res <- 0
  }
  res
}

# PAIRFST
# Pairwise FST among populations
# this functions computes pairwise FST between several populations and across several sites
# there is a single output FST value per comparison - the output is a triangular matrix
pairFST <- function(nPops, Pop_pi, Pool) {
  
  # create a matrix to store the FST values
  pairfst <- matrix(NA, ncol = nPops, nrow = nPops)
  
  ## Dealing with a matrix - single locus simulation
  ## ###############################################
  for(i in 1:(nPops-1)) {
    for(j in (i+1):nPops) {
      # apply the getFst() function to all the different pairwise combinations of the populations
      pairfst[i,j] <- getFst(Pop_pi[i,], Pop_pi[j,], Pool[i,], Pool[j,])
    }
  }
  
  # output FSt between pairs of populations
  pairfst
}

# POPSFST
# Pairwise FST among populations - across multiple loci
popsFST <- function(nPops, Pop_pi, Pool) {
  
  # check if the Pop_pi is in the correct (list) format
  if(class(Pop_pi) != "list")
    stop(paste("This function should be utilized with more than one locus. Please check"))
  
  ## Compute FST across a list - multiple locus simulation
  ## #####################################################
  listFST <- lapply(1:length(Pop_pi), function(locus) 
    pairFST(nPops = nPops, Pop_pi = Pop_pi[[locus]], Pool = Pool[["total"]][[locus]]))
  
  # output FST between pairs of populations
  listFST
}


# D - stat ----------------------------------------------------------------

# ABBA
# calculate the abba portion of the D-statistic
abba <- function(p1, p2, p3) {
  (p1 * (1 - p2) * (1 - p3)) + ((1 - p1) * p2 * p3)
}

# BABA
# calculate the baba portion of the D-statistic
baba <- function(p1, p2, p3) {
  ((1 - p1) * p2 * (1 - p3)) + (p1 * (1 - p2) * p3)
}


# D.STAT
# calculate D-statistic
D.stat <- function(ABBA, BABA) {
  (sum(ABBA) - sum(BABA)) / (sum(ABBA) + sum(BABA))
}

# Perform D-statistics analysis
# this functions calculates 3 different D-statistic values from pooled sequenced data -
# using the allelic frequencies for the ancestral allele 
D.statPool <- function(pop_pi) {
  
  ## Calculate D-stat with the crab population in site 1 as P3, wave pop in site 1 as P1 and wave pop in site 2 as P2
  ABBA <- abba(p1 = pop_pi[2,], p2 = pop_pi[4,], p3 = pop_pi[1,])
  BABA <- baba(p1 = pop_pi[2,], p2 = pop_pi[4,], p3 = pop_pi[1,])
  
  DStat1 <- D.stat(ABBA, BABA)
  
  ## Calculate D-stat with the crab population in site 1 as P3 and with the wave pop in site 1 as P1 
  # but with the crab pop at site 2 as P2
  ABBA <- abba(p1 = pop_pi[2,], p2 = pop_pi[3,], p3 = pop_pi[1,])
  BABA <- baba(p1 = pop_pi[2,], p2 = pop_pi[3,], p3 = pop_pi[1,])
  
  DStat2 <- D.stat(ABBA, BABA)
  
  ## Calculate D-stat with the wave population in site 2 as P3, wave pop in site 1 as P1 and crab pop in site 1 as P2
  ABBA <- abba(p1 = pop_pi[2,], p2 = pop_pi[1,], p3 = pop_pi[4,])
  BABA <- baba(p1 = pop_pi[2,], p2 = pop_pi[1,], p3 = pop_pi[4,])
  
  DStat3 <- D.stat(ABBA, BABA)
  
  ## Output the values of D-stat: With different populations as P3
  dstat <- c(DStat1, DStat2, DStat3)
  dstat
}


# master functions --------------------------------------------------------

# create the null output of the simulations
# when something fails, this will be the default output
# each possible sumstat will be set to NA
null.output <- function(nPops) {
  
  # when the models have four populations
  if(nPops == 4) {
    
    # create the output if the target number of loci was not reached
    out <- rep(list(NA), 64)
    
    # create a vector with the names of the output
    stats <- c("nPoly","nFilter","nLoci","Sf1","Sf2","Sf3","Sx1","Sx2","Sx3","Sx4","Sx5","SS1","SS2","SS3","Mean_Het1","Mean_Het2",
               "Mean_Het3","Mean_Het4","SD_Het1","SD_Het2","SD_Het3","SD_Het4","Mean_HetBet1","Mean_HetBet2","Mean_HetBet3",
               "Mean_HetBet4","Mean_HetBet5","Mean_HetBet6","SD_HetBet1","SD_HetBet2","SD_HetBet3","SD_HetBet4","SD_HetBet5",
               "SD_HetBet6","Mean_FST1","Mean_FST2","Mean_FST3","Mean_FST4","Mean_FST5","Mean_FST6","SD_FST1","SD_FST2","SD_FST3",
               "SD_FST4","SD_FST5","SD_FST6","FSTQ11","FSTQ12","FSTQ13","FSTQ14","FSTQ15","FSTQ16","FSTQ21","FSTQ22","FSTQ23",
               "FSTQ24","FSTQ25","FSTQ26","Dstat1","Dstat2","Dstat3","SD_dstat1","SD_dstat2","SD_dstat3")
    
    # add names to the output
    names(out) = stats
    
  } else { # for models with two populations
    
    # create the output if the target number of loci was not reached
    out <- rep(list(NA), 17)
    
    # create a vector with the names of the output
    stats <- c("nPoly", "nFilter", "nLoci", "Sf", "Sx1", "Sx2", "SS", "Mean_Het1", "Mean_Het2", "SD_Het1", "SD_Het2", 
               "Mean_HetBet", "SD_HetBet", "Mean_FST", "SD_FST", "FSTQ1", "FSTQ2")
    
    # add names to the output
    names(out) = stats
  }
  
  # output the result of the function - the null output of the models
  out
}

# ABC simulation of Pooled DNA sequencing
# this function utilizes all of the previous functions, except for the function that draws parameters from the prior
# it simulates pooled data and then computes various summary statistics from that data
poolABC <- function(parameters, model, nDip, nPops, size, nLoci, nSites, mutrate, mean, variance, 
                    minimum, maximum, min.minor = NA) {
  
  # create a variable to store the target number of loci
  target <- nLoci
  
  # since it is possible that, after removing sites (further ahead in this function), some loci are left without a single site
  # we should simulate more loci than required
  # create a variable to increase the number of simulated loci - call it extra
  extra <- target + 100
  
  # create the command line to run the scrm package
  # the command line varies according to the selected model
  if(model == "noMig") {
    # create the command line for a simple model with 2 populations and no migration
    commands <- cmdNoMig(parameters, nSites, nLoci = extra, nDip, nPops, mutrate)
    
  } else if (model == "simple2pops") {
    # create the command line for a simple model with 2 populations but with migration
    commands <- cmdSimple2pops(parameters, nSites, nLoci = extra, nDip, nPops, mutrate)
    
  } else if (model == "2pops") {
    # create the command line for a model with 2 populations with migration
    commands <- cmd2pops(parameters, nSites, nLoci = extra, nDip, nPops, mutrate)
    
  } else if (model == "mig2pops") {
    # create the command line for a model with 2 populations with migration
    commands <- mig2pops(parameters, nSites, nLoci = extra, nDip, nPops, mutrate)
    
  } else if (model == "star") {
    # create the command line for a star shaped model
    commands <- cmdStar(parameters, nSites, nLoci = extra, nDip, nPops, mutrate)
    
  } else if (model == "Parallel") {
    # create the command line for the parallel origin model
    commands <- cmdParallel(parameters, nSites, nLoci = extra, nDip, nPops, mutrate)
    
  } else if (model == "Single") {
    # create the command line for the single origin model
    commands <- cmdSingle(parameters, nSites, nLoci = extra, nDip, nPops, mutrate)
    
  } else {
    # if a correct model is not supplied as input for the function - stop and warn
    stop(paste("Model should be noMig, simple2pops, mig2pops, 2pops, star, Parallel or Single. Please check!"))
  }
  
  # run the scrm package and obtain the genotypes
  genotypes <- runSCRM(commands, nDip, nPops, model)
  
  # get the mean number of polymorphic sites
  nPoly <- mean(sapply(genotypes, ncol))
  
  # simulate number of reads
  reads <- simulateCoverage(mean, variance, genotypes)
  
  # remove sites with a depth of coverage above or below the defined threshold
  reads <- remove_by_reads(nLoci = extra, reads, minimum = minimum, maximum = maximum, genotypes = genotypes)
  
  # get the genotypes - without sites simulated with a coverage below or above the threshold
  genotypes <- lapply(reads, function(locus) locus[[2]])
  
  # get the reads - without sites simulated with a coverage below or above the threshold
  reads <- lapply(reads, function(locus) locus[[1]])
  
  # ensure that each entry is a matrix
  reads <- lapply(reads, function(locus) matrix(locus, nrow = nPops))
  
  # simulate individual contribution to the total number of reads
  indContribution <- lapply(1:extra, function(locus)
    popsReads(list_np = size, coverage = reads[[locus]], pError = parameters["Pool_Error"]))
  
  # simulate the number of Ancestral reads
  ancestral <- lapply(1:extra, function(locus)
    numberAncestralPop(genotypes = genotypes[[locus]], indContribution = indContribution[[locus]], 
                       size = size, error = parameters["Error"]))
  
  # simulate pooled sequencing data
  pool <- poolPops(nPops = nPops, nLoci = extra, indContribution = indContribution, readsAncestral = ancestral)
  
  # use an lapply to ensure that the ancestral allele of the simulations is also the major allele - do this for each locus
  pool <- lapply(1:extra, function(locus) 
    majorSims(ancestral = pool[["Ancestral"]][[locus]], derived = pool[["Derived"]][[locus]], 
              coverage = pool[["Total"]][[locus]], min.minor = 2))
  
  # convert the pool list back to the previous format 
  # one entry for ancestral matrices, one for derived and a final one for total matrices
  pool <- list(ancestral = lapply(pool, function(locus) locus[["Ancestral"]]), 
               derived = lapply(pool, function(locus) locus[["Derived"]]), 
               total = lapply(X = pool, function(locus) locus[["Total"]]))
  
  # randomly select "nLoci" entries from the lists containing the pooled data
  pool <- pickLoci(commands, nPops, nLoci = nLoci, extra, major = pool$ancestral, minor = pool$derived, coverage = pool$total)
  
  # it is possible that some loci were eliminated during the previous steps
  # thus, it might be necessary to update the nLoci variable
  nLoci <- length(pool[["ancestral"]])
  
  # when the simulation did not produce the required number of loci
  if(nLoci < target) {
    
    # set all the values of the possible sumstats to NA
    output <- null.output(nPops = nPops)
    
  } else {
    
    # compute the fraction of sites showing a fixed difference between the populations
    Sf <- lapply(1:nLoci, function(locus) 
      fixed(ancestral = pool[["ancestral"]][[locus]], total = pool[["total"]][[locus]], nPops))
    # combine the previous list into a single matrix - where each column is a different comparison
    Sf <- do.call(rbind, Sf)
    # get the mean, across all loci, for the fraction of sites showing a fixed difference
    Sf <- colMeans(Sf, na.rm = TRUE)
    
    # calculate the fraction of sites showing an exclusive polymorphism to a given population
    Sx <- lapply(1:nLoci, function(locus) 
      exclusive(ancestral = pool[["ancestral"]][[locus]], total = pool[["total"]][[locus]], nPops))
    # combine the previous list into a single matrix - where each column is a different comparison
    Sx <- do.call(rbind, Sx)
    # compute the mean (across loci) for each population
    Sx <- colMeans(Sx, na.rm = TRUE)
    
    # compute the fraction of sites with a polymorphism shared between the populations
    SS <- lapply(1:nLoci, function(locus) 
      shared(ancestral = pool[["ancestral"]][[locus]], total = pool[["total"]][[locus]], nPops))
    # combine the previous list into a single matrix - where each column is a different comparison
    SS <- do.call(rbind, SS)
    # get the mean, across all loci, for the fraction of sites shared polymorphism between the populations
    SS <- colMeans(SS, na.rm = TRUE)
    
    # Calculate population allelic frequency
    Pop_Pi <- calculatePi(listPool = pool, nLoci = nLoci)
    pool <- Pop_Pi[[2]]
    Pop_Pi <- Pop_Pi[[1]]
    
    # get the mean number of polymorphic sites - after filtering
    nFilter <- mean(sapply(Pop_Pi, ncol))
    
    # Compute the mean expected heterozygosity for each population and locus
    ExpHet <- meanExpected_Het(Pop_Pi)
    # Compute the mean (across loci) for each population
    PopHet <- do.call(rbind,ExpHet)
    # Compute mean expected heterozygosity for each population across all loci
    HetMean <- colMeans(PopHet, na.rm = TRUE)
    # And the standard deviation
    HetSD <- apply(PopHet, 2, sd, na.rm = TRUE)
    
    # Compute the heterozygosity between populations
    HetBetween <- Het_Between(Pop_Pi)
    # Compute the mean (across loci) for each pairwise comparison
    PopBetHet <- do.call(rbind,HetBetween)
    MeanHetBet <- colMeans(PopBetHet, na.rm = TRUE)
    # And the standard deviation
    SDHetBet <- apply(PopBetHet, 2, sd, na.rm = TRUE)
    
    # Calculate the FST value between populations
    FST <- popsFST(nPops, Pop_Pi, pool)
    # The previous returns a matrix where each pairwise comparison has a FST value and the rest is NAs
    # This matrix can be reduced to a vector
    FST <- lapply(FST, FUN = function(x) {
      x[!is.na(x)]})
    # combine all entries into a single matrix
    FST <- do.call(rbind,FST)
    # replace negative values by zero
    FST[which(FST < 0)] = 0
    # Then, compute the mean FST value between pops
    MeanFST <- colMeans(FST, na.rm = TRUE)
    # And the standard deviation
    SDFST <- apply(FST, MARGIN = 2, sd)
    # calculate the 5% and the 95% quantiles for the FST distribution
    FSTQ1 <- apply(FST, MARGIN = 2, function(col) unname(quantile(col, probs = 0.05)))
    FSTQ2 <- apply(FST, MARGIN = 2, function(col) unname(quantile(col, probs = 0.95)))
    
    # compute Dstat values - if 4 populations were used
    if (nPops == 4) {
      # calculate DSTAT values over a list - each entry is a different locus
      dstat <- lapply(Pop_Pi, function(pi) D.statPool(pi))
      # combine all the values into a single matrix
      tempdstat <- do.call(rbind,dstat)
      # compute the mean value for a single simulation
      dstat <- colMeans(tempdstat, na.rm = TRUE)
      # And the standard deviation across locus
      SD_dstat <- apply(tempdstat, MARGIN = 2 , sd, na.rm = TRUE)
      
      # create the output
      output <- list(nPoly, nFilter, nLoci, Sf, Sx, SS, HetMean, HetSD, MeanHetBet, SDHetBet, MeanFST, SDFST, FSTQ1, FSTQ2, 
                     dstat, SD_dstat)
      names(output) = c("nPoly", "nFilter", "nLoci", "Sf", "Sx", "SS", "Mean_Het", "SD_Het", "Mean_HetBet", "SD_HetBet", 
                        "Mean_FST", "SD_FST", "FSTQ1", "FSTQ2", "Dstat", "SD_dstat")
      
    } else { # if we are working with models that don't have 4 populations
      
      # create the output
      output <- list(nPoly, nFilter, nLoci, Sf, Sx, SS, HetMean, HetSD, MeanHetBet, SDHetBet, MeanFST, SDFST, FSTQ1, FSTQ2)
      names(output) = c("nPoly", "nFilter", "nLoci", "Sf", "Sx", "SS", "Mean_Het", "SD_Het", "Mean_HetBet", "SD_HetBet", 
                        "Mean_FST", "SD_FST", "FSTQ1", "FSTQ2")
    }
  }
  
  # output the final result of the function
  output
}

# Perform ABC simulation with pooled sequence data
# this is a master function - it draws parameters from the priors, simulates pooled sequence data and computes sumstats
# the output consists of both the drawn parameters and the calculated summary statistics 
poolStats <- function(model, nDip, nPops, size, nLoci, nSites, mutrate, mean, variance, minimum, maximum, min.minor = NA, NE, 
                      ratio, split, pool, seq, CW = NA, WC = NA, Mb = NA, Md = NA, site = NA, anc = NA, bT = NA, bCW = NA, 
                      bWC = NA, pMd = NA) {
  
  # check if the input is correct - when using the noMig, mig2pops, simple2pops or 2pops models, the nPops input should be 2
  if(model %in% c("noMig", "simple2pops", "mig2pops", "2pops") &  nPops != 2)
    stop(paste("The selected model should contain only two populations. Please check"))
  
  # check if the input is correct - the list "size" should contain one different entry per population
  if(class(size) == "list" & length(size) != nPops)
    stop(paste("The list input - size - should have one entry per population. Please check"))
  
  # check if the input is correct - the vectors "mean" and "variance" should contain one different entry per population
  if(length(mean) != nPops | length(variance) != nPops)
    stop(paste("The mean and variances inputs should have one entry per population. Please check"))
  
  # create a dataframe with the parameters 
  parameters <- createParams(NE, ratio, split, pool, seq, CW, WC, Mb, Md, site, anc, bT, bCW, bWC, pMd, model)
  
  # use the parameters dataframe as input for the poolABC function
  SumStats <- poolABC(parameters, model, nDip, nPops, size, nLoci, nSites, mutrate, mean, variance, minimum, maximum, min.minor)
  
  # output the summary statistics and the parameters that led to those summary statistics
  c(parameters, SumStats)
}
