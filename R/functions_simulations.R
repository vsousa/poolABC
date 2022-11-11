#' Draw parameters from the priors
#'
#' This function creates a named vector of parameters that can be used as input
#' in the command line of the scrm package. Please note that this function needs
#' to be adjusted if you wish to test the effect of different prior
#' distributions.
#'
#' @param Nref The minimum and maximum value of the uniform distribution for the
#'   effective population size of the reference population (Nref).
#' @param ratio The minimum and maximum value of the distribution from which the
#'   relative size of the present-day and ancestral populations are drawn. The
#'   size of these populations is set as a ratio of the size of the Nref
#'   population. All of these ratios are drawn from a log10 uniform
#'   distribution.
#' @param split The minimum and maximum values, at the 4Nref scale, of the
#'   uniform distribution from which the values of the times of the split events
#'   are draw. Both the time of the recent split event and the distance between
#'   the two split events are drawn from this distribution.
#' @param pool The minimum and maximum values of the uniform distribution from
#'   which the value of the error associated with DNA pooling is drawn. More
#'   specifically, this value is related with the unequal individual
#'   contribution to the pool.
#' @param seq The minimum and maximum values of the uniform distribution from
#'   which the value of the error associated with DNA sequencing is drawn. This
#'   parameter should be supplied as a decimal number between zero and one.
#' @param CW The minimum and maximum value of the uniform distribution from
#'   which the migration rate between the two divergent ecotypes inhabiting the
#'   same location is drawn. We consider that this parameter is drawn on a m
#'   scale. This is the migration rate from ecotype C to ecotype W.
#' @param WC The minimum and maximum value of the uniform distribution from
#'   which the migration rate between the two divergent ecotypes inhabiting the
#'   same location is drawn. We consider that this parameter is drawn on a m
#'   scale. This is the migration rate from ecotype W to ecotype C.
#' @param CC The minimum and maximum value of the uniform distribution from
#'   which the migration rate between similar ecotypes inhabiting different
#'   locations is drawn. We consider that this parameter is drawn on a m scale.
#'   This is the migration between the two C ecotypes at two different
#'   locations.
#' @param WW The minimum and maximum value of the uniform distribution from
#'   which the migration rate between similar ecotypes inhabiting different
#'   locations is drawn. We consider that this parameter is drawn on a m scale.
#'   This is the migration between the two W ecotypes at two different
#'   locations.
#' @param ANC The minimum and maximum value of the uniform distribution from
#'   which the migration rate between the two ancestral populations is drawn. We
#'   consider that this parameter is drawn on a m scale.
#' @param bT The minimum and maximum values of the distribution from which the
#'   proportion of the simulated loci where no migration occurs between
#'   divergent ecotypes is drawn. The maximum value should not be higher than
#'   one.
#' @param bCW The minimum and maximum values of the distribution from which the
#'   proportion of the simulated loci where no migration occurs from the C
#'   ecotype towards the W ecotype is drawn. The maximum value should not be
#'   higher than one.
#' @param bWC The minimum and maximum values of the distribution from which the
#'   proportion of the simulated loci where no migration occurs from the W
#'   ecotype towards the C ecotype is drawn. The maximum value should not be
#'   higher than one.
#' @param model Either "2pops", "Single" or "Parallel" indicating for which
#'   model should parameters be drawn.
#' @param digits An optional integer indicating the number of decimal places to
#'   use when rounding certain parameters. The default is five.
#'
#' @return a vector with one named entry per relevant parameter. Each entry is
#'   the sampled value from the prior for that particular parameter.
#'
#' @examples
#' # for a model with two populations
#' createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250), seq = c(0.0001, 0.001),
#' split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3), bT = c(0, 0.2), model = "2pops")
#'
#' # for a single origin scenario
#' createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250), seq = c(0.0001, 0.001),
#' split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3), CC =  c(1e-13, 1e-3),
#' WW = c(1e-13, 1e-3), ANC = c(1e-13, 1e-3), bT = c(0, 0.2), bCW = c(0, 0.5),
#' bWC = c(0, 0.5), model = "Single")
#'
#' @export
createParams <- function(Nref, ratio, split, pool, seq, CW, WC, CC = NA, WW = NA, ANC = NA, bT, bCW = NA, bWC = NA,
                         model, digits = 5) {

  # check if the input is correct - the selected model should be one of the following
  if(model %in% c("2pops", "Single", "Parallel") == FALSE)
    stop("The selected model should be either 2pops, Single or Parallel. Please check")

  # draw a value for the population size of the most ancestral population - it's also the Nref
  Nrf <- stats::runif(n = 1, min = Nref[1], max = Nref[2])
  # Draw the values for the extant and ancestral population sizes
  # Values are drawn as a ratio of the Ne size - a value of 2 means that the population is twice the size of the Ne
  N1 <- round(stats::runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  N2 <- round(stats::runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  # set the values in the natural scale
  N1 <- 10^N1; N2 <- 10^N2

  # Draw the split times from a uniform distribution - times are drawn on a 4Nref scale
  Split <- stats::runif(n = 1, min = split[1], max = split[2])
  # Draw the errors for the pooling parameter
  Pool_Error <- stats::runif(n = 1, min = pool[1], max = pool[2])
  # and the sequencing error
  Error <- stats::runif(n = 1, min = seq[1], max = seq[2])

  # draw the migration rate (at the m scale) - this corresponds to the migration between different ecotypes in the same site
  # the migration rate from crab to wave is drawn as mCW - for the first site
  mCW1 <- stats::runif(n = 1, min = CW[1], max = CW[2])
  # and the migration rate from wave to crab is drawn as mWC - for the first site
  mWC1 <- stats::runif(n = 1, min = WC[1], max = WC[2])

  # proportion of the genome without migration - total barrier
  total <- stats::rbeta(n = 1, shape1 = 1, shape2 = 10)
  # replace values below the minimum threshold with the minimum
  total[total < 1e-2] <- bT[1]
  # and values above the maximum threshold with the maximum
  total[total > bT[2]] <- bT[2]

  # stop the function if we are working with the two-population model
  if(model == "2pops") {
    # assume that the proportion of the genome with unrestricted migration is 1 - proportion without migration
    pMig <- 1 - total
    # create the parameters vector for this particular model
    parameters <- c(Nrf, N1, N2, Split, Pool_Error, Error, mCW1, mWC1, pMig, total)
    # add names to the entries of the vector
    names(parameters) <- c("Nref", "N1", "N2", "Split", "PoolError", "SeqError", "mCW", "mWC", "pM", "pNO")
    # stop the function and output the parameters vector
    stop(return(parameters))
  }

  # draw the values for the remaining present-day population sizes
  N3 <- round(stats::runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  N4 <- round(stats::runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  # set the values in the natural scale
  N3 <- 10^N3; N4 <- 10^N4

  # draw the values for the ancestral population sizes
  # values are drawn as a ratio of the Ne size - a value of 2 means that the population is twice the size of the Ne
  NA1 <- round(stats::runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  NA2 <- round(stats::runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  # set the values in the natural scale
  NA1 <- 10^NA1; NA2 <- 10^NA2

  # draw the second split time from a uniform distribution - times are drawn on a 4Ne scale
  Dsplit <- stats::runif(n = 1, min = split[1], max = split[2])

  # draw the migration rate (at the m scale) - this corresponds to the migration between different ecotypes in the same site
  # the migration rate from crab to wave is drawn as mCW - for the second site
  mCW2 <- stats::runif(n = 1, min = CW[1], max = CW[2])
  # and the migration rate from wave to crab is drawn as mWC - for the second site
  mWC2 <- stats::runif(n = 1, min = WC[1], max = WC[2])

  # if required, draw additional migration rates
  if(!any(is.na(CC)))
    mCC <- stats::runif(n = 1, min = CC[1], max = CC[2]) # between crab ecotypes in different locations
  else
    mCC <- NA  # set mCC to NA

  # if required, draw additional migration rates
  if(!any(is.na(WW)))
    mWW <- stats::runif(n = 1, min = WW[1], max = WW[2]) # between wave ecotypes in different locations
  else
    mWW <- NA  # set mWW to NA

  # if required, draw additional migration rates
  if(!any(is.na(ANC)))
    mAA <- stats::runif(n = 1, min = ANC[1], max = ANC[2])  # migration rates between the ancestral populations
  else
    mAA <- NA # set mAA to NA

  # if required, draw the proportion of the genome without migration
  if(!any(is.na(bCW))) {
    # from the crab to the wave ecotype
    pCW <- stats::rbeta(n = 1, shape1 = 1, shape2 = 10)
    # replace values below the minimum threshold with the minimum
    pCW[pCW < 1e-2] <- bCW[1]
    # and values above the maximum threshold with the maximum
    pCW[pCW > bCW[2]] <- bCW[2]

  } else {

    # set pCW to NA
    pCW <- NA
  }

  # if required, draw the proportion of the genome without migration
  if(!any(is.na(bWC))) {
    # from the wave to the crab ecotype
    pWC <- stats::rbeta(n = 1, shape1 = 1, shape2 = 10)
    # replace values below the minimum threshold with the minimum
    pWC[pWC < 1e-2] <- bWC[1]
    # and values above the maximum threshold with the maximum
    pWC[pWC > bWC[2]] <- bWC[2]

  } else {

    # set pWC to NA
    pWC <- NA
  }

  # assume that the proportion of the genome with unrestricted migration is 1 - proportion without migration
  pMig <- 1 - sum(c(total, pCW, pWC), na.rm = TRUE)

  # create the parameters vector for the single and parallel models
  parameters <- c(Nrf, N1, N2, N3, N4, NA1, NA2, Split, Dsplit, Pool_Error, Error, mCW1, mCW2, mWC1, mWC2, mCC, mWW, mAA,
                  pMig, pCW, pWC, total)

  # add names to the entries of the vector
  names(parameters) <- c("Nref", "N1", "N2", "N3", "N4", "NA1", "NA2", "Split", "Dsplit", "PoolError", "SeqError", "mCW1", "mCW2",
                         "mWC1", "mWC2", "mCC", "mWW", "mAA", "pM", "pCW", "pWC", "pNO")

  # remove any parameters that are set as NA
  parameters <- parameters[!is.na(parameters)]

  # output the parameters vector
  parameters
}


#' Create SCRM command line for a model with two populations
#'
#' This function creates a command line tailored for an isolation with migration
#' model with two populations. The command line can then be fed to the scrm
#' package to run the model.
#'
#' @param parameters A vector where each entry corresponds to a different
#'   parameter, e.g. one entry is the size of the reference population, another
#'   is the time of recent split, etc. Please note that this functions depends
#'   on the ordering of the parameters in the vector and thus, it should only be
#'   used with a vector created with the `createParams` function.
#' @param nSites An integer representing the number of base pairs that each
#'   locus should have.
#' @param nLoci An integer that represents how many independent loci should be
#'   simulated.
#' @param nDip An integer representing the total number of diploid individuals
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this. Also note that this is the total
#'   number of diploid individuals and this function will distribute the
#'   individuals equally by the two populations.
#' @param mutrate A number representing the mutation rate assumed for the
#'   simulations.
#' @param extra is a logical value indicating whether the required number of
#'   loci should be enforced. The default is FALSE but, if set to TRUE, then
#'   additional loci will be simulated. These additional loci are simulated to
#'   try to have sufficient loci to keep the required number of loci after
#'   filtering.
#'
#' @return a character vector with two entries. The first entry is the scrm
#'   command line for the loci without any barriers against migration, while the
#'   second entry is the scrm command line for the loci without migration
#'   between divergent ecotypes.
#'
#' @examples
#' # create a vector with parameter values for a two populations model
#' params <- createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250),
#' seq = c(0.0001, 0.001), split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3),
#' bT = c(0, 0.2), model = "2pops")
#'
#' # create the command line for the scrm package
#' cmd2pops(parameters = params, nSites = 2000, nLoci = 100, nDip = 100, mutrate = 2e-8)
#'
#' @export
cmd2pops <- function(parameters, nSites, nLoci, nDip, mutrate, extra = FALSE) {

  # this function is intended to be used with a two-population model
  nPops <- 2

  # Read the vector with the parameters and assign each parameter to the correct command name
  Ne <- parameters[1]
  # set the relative size of each population
  N1 <- parameters[2]; N2 <- parameters[3]

  # get the proportion of loci with migration - no barriers against migration between the different ecotypes
  pM <- parameters[9]
  # and the proportion of loci without any migration - total barrier against migration between the different ecotypes
  pNO <- parameters[10]

  # get the migration rates - between different ecotypes at the same site - this is the m value on the M = 4N0m formula
  # from the crab ecotype to the wave ecotype
  mCW <-  parameters[7]
  # from the wave ecotype to the crab ecotype
  mWC <-  parameters[8]

  # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
  # and REMEMBER that M = 4N0m
  # set the migration rates to the scale of Nref - between different ecotypes at the same site
  mig_CW <- 4*Ne*mCW # -m 2 1 mig_CW

  # the migration from wave to crab is parametrized as a ratio of the migration from crab to wave
  # so we need to multiply the migration from crab to wave by this ratio
  mig_WC <- 4*Ne*mWC # -m 1 2 mig_WC

  # get the time of the split event
  split <- round(parameters[4], digits = 3)

  # Compute the value of theta
  mutrate_locus <- nSites*mutrate
  theta <- 4*Ne*mutrate_locus
  # Create a vector with information about how many haplotypes are sampled from each population
  n <- c(rep(nDip/nPops, times = nPops))

  # use a multinomial distribution to get the number of loci simulated under each category
  lociTotal <- as.vector(stats::rmultinom(n = 1, size = nLoci, c(pM, pNO)))

  # if extra is TRUE, then more loci than required per category will be simulated
  if(extra == TRUE) {
    # save the required number of simulated loci per category
    targetLoci <- lociTotal
    # simulate more 25 loci per category
    lociTotal <- lociTotal + 25
  }

  # cheat code: pop1 - crab in site 1; pop2 - wave in site 1; pop3 - crab in site 2; pop4 - wave in site 2
  # create command line with no barriers to migration
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and
  # nrep: the number of independent loci that will be produced
  with.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[1], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                    # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                    "-n 1", N1, "-n 2", N2,
                    # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
                    # set the migration rate between different ecotypes inhabiting the same location
                    "-m 2 1", mig_CW, "-m 1 2", mig_WC,
                    # now, set the migration rate right before the split event to zero by using the switch:
                    # -eM <t> <M>: assume a symmetric migration rate of M/(npop-1) at time t.
                    "-eM", split, "0",
                    # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                    # finally, set the size of the ancestral pop equal to the size of the reference population with:
                    # -eN <t> <n>: set the size of all populations to n*N0 at time t.
                    "-ej", split, "2 1 -eN", split, 1)

  # create a command line for the loci without any migration (between the different ecotypes)
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and
  # nrep: the number of independent loci that will be produced
  without.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[2], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "), "0",
                       # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                       "-n 1", N1, "-n 2", N2,
                       # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                       # finally, set the size of the ancestral pop equal to the size of the reference population with:
                       # -eN <t> <n>: set the size of all populations to n*N0 at time t.
                       "-ej", split, "2 1 -eN", split, 1)


  # combine the two different types of commands
  cmd_2pops <- c(with.mig, without.mig)

  # if extra is equal to TRUE, then we simulated more loci than required
  if(extra == TRUE)
    # include the required number of loci per category in the output
    cmd_2pops <- list(commands = cmd_2pops, targetLoci = targetLoci)

  # output the command line for the two population model
  cmd_2pops
}

#' Create SCRM command line for a single origin scenario
#'
#' This function creates a command line tailored for a scenario of single origin
#' to explain ecotype formation. The command line can then be fed to the scrm
#' package to run the model.
#'
#' For convenience, imagine we have two divergent ecotypes, named C and W. This
#' model assumes that the first population corresponds to the C ecotype at the
#' first location, the second population to the C ecotype in the first location,
#' the third population to the W ecotype in the second location and the fourth
#' population to the W ecotype in the second location.
#'
#' @param parameters A vector where each entry corresponds to a different
#'   parameter, e.g. one entry is the size of the reference population, another
#'   is the time of recent split, etc. Please note that this functions depends
#'   on the ordering of the parameters in the vector and thus, it should only be
#'   used with a vector created with the `createParams` function.
#' @param nSites An integer representing the number of base pairs that each
#'   locus should have.
#' @param nLoci An integer that represents how many independent loci should be
#'   simulated.
#' @param nDip An integer representing the total number of diploid individuals
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this. Also note that this is the total
#'   number of diploid individuals and this function will distribute the
#'   individuals equally by the two populations.
#' @param mutrate A number representing the mutation rate assumed for the
#'   simulations.
#' @param extra is a logical value indicating whether the required number of
#'   loci should be enforced. The default is FALSE but, if set to TRUE, then
#'   additional loci will be simulated. These additional loci are simulated to
#'   try to have sufficient loci to keep the required number of loci after
#'   filtering.
#'
#' @return a character vector with four entries. The first entry is the scrm
#'   command line for the loci without any barriers against migration. The
#'   second entry is the command line for the loci without migration from the C
#'   towards the W ecotype. The third entry is command line for the loci without
#'   migration from the W towards the C ecotype and the last entry is the scrm
#'   command line for the loci without migration between divergent ecotypes.
#'
#' @examples
#' # create a vector with parameter values for the single origin scenario
#' params <- createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250),
#' seq = c(0.0001, 0.001), split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3),
#' CC =  c(1e-13, 1e-3), WW = c(1e-13, 1e-3), ANC = c(1e-13, 1e-3), bT = c(0, 0.2),
#' bCW = c(0, 0.5), bWC = c(0, 0.5), model = "Single")
#'
#' # create the command line for the scrm package
#' cmdSingle(parameters = params, nSites = 2000, nLoci = 100, nDip = 400, mutrate = 2-8)
#'
#' @export
cmdSingle <- function(parameters, nSites, nLoci, nDip, mutrate, extra = FALSE) {

  # this function is intended to be used with a four-population model
  nPops <- 4

  # read the vector with the parameters and assign each parameter to the correct variable name
  Ne <- parameters[1]
  # set the relative size of each population - for the extant populations
  N1 <- parameters[2]; N2 <- parameters[3]; N3 <- parameters[4]; N4 <- parameters[5]
  # and the ancient populations
  NA1 <- parameters[6]; NA2 <- parameters[7]

  # get the proportion of loci with migration - no barriers against migration between the different ecotypes
  pM <- parameters[19]
  # get the proportion of loci without migration - from the crab to the wave ecotype at the same location
  pCW <- parameters[20]
  # get the proportion of loci without migration - from the wave to the crab ecotype at the same location
  pWC <- parameters[21]
  # and the proportion of loci without any migration - total barrier against migration between the different ecotypes
  pNO <- parameters[22]

  # get the migration rates - between different ecotypes at the same site - this is the m value on the M = 4N0m formula
  # from the crab ecotype to the wave ecotype - at the first site
  mCW1 <- parameters[12]
  # from the crab ecotype to the wave ecotype - at the second site
  mCW2 <- parameters[13]
  # from the wave ecotype to the crab ecotype - at the fist site
  mWC1 <- parameters[14]
  # from the wave ecotype to the crab ecotype - at the second site
  mWC2 <- parameters[15]

  # between crab populations inhabiting different locations
  mCC <- parameters[16]
  # between wave populations inhabiting different locations
  mWW <- parameters[17]
  # and between the two ancestral populations
  mAA <- parameters[18]

  # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
  # REMEMBER that M = 4N0m
  # set the migration rates to the scale of Nref - between different ecotypes at the same site
  # from the crab ecotype to the wave ecotype - at the first site
  mig_CW1 <- 4*Ne*mCW1 # -m 3 1 mig_CW
  # from the crab ecotype to the wave ecotype - at the second site
  mig_CW2 <- 4*Ne*mCW2 # -m 4 2 mig_CW
  # from the wave ecotype to the crab ecotype - at the first site
  mig_WC1 <- 4*Ne*mWC1 # -m 1 3 mig_WC
  # from the wave ecotype to the crab ecotype - at the second site
  mig_WC2 <- 4*Ne*mWC2 # -m 2 4 mig_WC

  # between crab populations at different locations
  mig_CC <- 4*Ne*mCC
  # between wave populations at different locations
  mig_WW <- 4*Ne*mWW
  # and between the two ancestral populations
  mig_AA <- 4*Ne*mAA

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

  # create variables to ensure that the changes in migration rates occur at different times than the split time
  # a variable to inform when does migration start between ancestral populations
  tmAA <- Rsplit + 0.0001
  # create a variable to inform when does migration stop before the recent split
  # if Rsplit is zero, we can not subtract something from it
  if(Rsplit != 0) {
    # when Rsplit is not zero, set the end of the migration before the split
    tmRS <- Rsplit - 0.0001
  } else {
    # when Rsplit is zero, set the end of the migration at the split
    tmRS <- Rsplit
  }

  # create also a variable to inform when does migration start between ancestral populations
  tmAA <- Rsplit + 0.0001

  # use a multinomial distribution to get the number of loci simulated under each category
  lociTotal <- as.vector(stats::rmultinom(n = 1, size = nLoci, c(pM, pCW, pWC, pNO)))

  # if extra is TRUE, then more loci than required per category will be simulated
  if(extra == TRUE) {
    # save the required number of simulated loci per category
    targetLoci <- lociTotal
    # simulate more 25 loci per category
    lociTotal <- lociTotal + 25
  }

  # cheat code: pop1 - crab in site 1; pop2 - crab in site 2; pop3 - wave in site 1; pop4 - wave in site 2

  # create command line with no barriers to migration
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and
  # nrep: the number of independent loci that will be produced
  with.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[1], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                    # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                    "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                    # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
                    # set the migration rate between different ecotypes inhabiting the same location
                    "-m 3 1", mig_CW1, "-m 4 2", mig_CW2, "-m 1 3", mig_WC1, "-m 2 4", mig_WC2,
                    # set the migration rate between the same ecotypes inhabiting different locations
                    "-m 2 1", mig_CC, "-m 1 2", mig_CC, "-m 4 3", mig_WW, "-m 3 4", mig_WW,
                    # set the migration, between all populations, to zero - immediately before the recent split
                    "-eM", tmRS, "0",
                    # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                    # looking backwards in time, it moves all lines from population j into population i at time t
                    # Migration rates into population j are set to 0 for the time further back into the past
                    "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                    # set the size of the ancient populations
                    # -en <t> <i> <n>: Set the size of population i to n*N0 at time t.
                    "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                    # set the migration rate between the ancestral populations
                    "-em", tmAA, "2 3", mig_AA, "-em", tmAA, "3 2", mig_AA,
                    # set the migration, between all populations, to zero - immediately at the ancient split
                    "-eM", Asplit, "0",
                    # add a split event - this event creates the two ancestral populations
                    "-ej", Asplit, "2 3",
                    # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                    # -eN <t> <n>: set the size of all populations to n*N0 at time t.
                    "-eN", Asplit, 1)

  # create command line with no migration from the crab to the wave ecotype within the same location
  no.mCW <- paste(paste(nDip*2, collapse = " "), lociTotal[2], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                  # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0
                  "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                  # set the migration rate from wave to crab (looking forward in time) inhabiting the same location
                  "-m 3 1 0 -m 4 2 0 -m 1 3", mig_WC1, "-m 2 4", mig_WC2,
                  # set the migration rate between the same ecotypes inhabiting different locations
                  "-m 2 1", mig_CC, "-m 1 2", mig_CC, "-m 4 3", mig_WW, "-m 3 4", mig_WW,
                  # set the migration, between all populations, to zero - immediately before the recent split
                  "-eM", tmRS, "0",
                  # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                  "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                  # set the size of the ancient populations
                  "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                  # set the migration rate between the ancestral populations
                  "-em", tmAA, "2 3", mig_AA,
                  # set the migration, between all populations, to zero - immediately at the ancient split
                  "-eM", Asplit, "0",
                  # add a split event - this event creates the two ancestral populations
                  "-ej", Asplit, "2 3",
                  # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                  "-eN", Asplit, 1)

  # create command line with no migration from the wave to the crab ecotype within the same location
  no.mWC <- paste(paste(nDip*2, collapse = " "), lociTotal[3], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                  # set the size of the present day populations
                  "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                  # set the migration rate from crab to wave (looking forward in time) inhabiting the same location
                  "-m 3 1", mig_CW1, "-m 4 2", mig_CW1, "-m 1 3 0 -m 2 4 0",
                  # set the migration rate between the same ecotypes inhabiting different locations
                  "-m 2 1", mig_CC, "-m 1 2", mig_CC, "-m 4 3", mig_WW, "-m 3 4", mig_WW,
                  # set the migration, between all populations, to zero - immediately before the recent split
                  "-eM", tmRS, "0",
                  # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                  "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                  # set the size of the ancient populations
                  "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                  # set the migration rate between the ancestral populations
                  "-em", tmAA, "3 2", mig_AA,
                  # set the migration, between all populations, to zero - immediately at the ancient split
                  "-eM", Asplit, "0",
                  # add a split event - this event creates the two ancestral populations
                  "-ej", Asplit, "2 3",
                  # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                  "-eN", Asplit, 1)

  # create a command line for the loci without any migration (between the different ecotypes)
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and
  # nrep: the number of independent loci that will be produced
  without.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[4], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "), "0",
                       # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                       "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                       # set the migration rate between the same ecotypes inhabiting different locations
                       "-m 2 1", mig_CC, "-m 1 2", mig_CC, "-m 4 3", mig_WW, "-m 3 4", mig_WW,
                       # set the migration, between all populations, to zero - immediately before the recent split
                       "-eM", tmRS, "0",
                       # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                       "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                       # set the size of the ancient populations
                       "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                       # set the migration, between all populations, to zero - immediately at the ancient split
                       "-eM", Asplit, "0",
                       # add a split event - this event creates the two ancestral populations
                       "-ej", Asplit, "2 3",
                       # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                       "-eN", Asplit, 1)

  # combine the two different types of commands
  cmdSingle <- c(with.mig, no.mCW, no.mWC, without.mig)

  # if extra is equal to TRUE, then we simulated more loci than required
  if(extra == TRUE)
    # include the required number of loci per category in the output
    cmdSingle <- list(commands = cmdSingle, targetLoci = targetLoci)

  # output the command line for the single origin model
  cmdSingle
}


#' Create SCRM command line for a parallel origin scenario
#'
#' This function creates a command line tailored for a scenario of parallel
#' origin to explain ecotype formation. The command line can then be fed to the
#' scrm package to run the model.
#'
#' For convenience, imagine we have two divergent ecotypes, named C and W. This
#' model assumes that the first population corresponds to the C ecotype at the
#' first location, the second population to the W ecotype in the first location,
#' the third population to the C ecotype in the second location and the fourth
#' population to the W ecotype in the second location.
#'
#' @param parameters A vector where each entry corresponds to a different
#'   parameter, e.g. one entry is the size of the reference population, another
#'   is the time of recent split, etc. Please note that this functions depends
#'   on the ordering of the parameters in the vector and thus, it should only be
#'   used with a vector created with the `createParams` function.
#' @param nSites An integer representing the number of base pairs that each
#'   locus should have.
#' @param nLoci An integer that represents how many independent loci should be
#'   simulated.
#' @param nDip An integer representing the total number of diploid individuals
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this. Also note that this is the total
#'   number of diploid individuals and this function will distribute the
#'   individuals equally by the two populations.
#' @param mutrate A number representing the mutation rate assumed for the
#'   simulations.
#' @param extra is a logical value indicating whether the required number of
#'   loci should be enforced. The default is FALSE but, if set to TRUE, then
#'   additional loci will be simulated. These additional loci are simulated to
#'   try to have sufficient loci to keep the required number of loci after
#'   filtering.
#'
#' @return a character vector with four entries. The first entry is the scrm
#'   command line for the loci without any barriers against migration. The
#'   second entry is the command line for the loci without migration from the C
#'   towards the W ecotype. The third entry is command line for the loci without
#'   migration from the W towards the C ecotype and the last entry is the scrm
#'   command line for the loci without migration between divergent ecotypes.
#'
#' @examples
#' # create a vector with parameter values for the parallel origin scenario
#' params <- createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250),
#' seq = c(0.0001, 0.001), split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3),
#' CC =  c(1e-13, 1e-3), WW = c(1e-13, 1e-3), ANC = c(1e-13, 1e-3), bT = c(0, 0.2),
#' bCW = c(0, 0.5), bWC = c(0, 0.5), model = "Parallel")
#'
#' # create the command line for the scrm package
#' cmdParallel(parameters = params, nSites = 2000, nLoci = 100, nDip = 400, mutrate = 2-8)
#'
#' @export
cmdParallel <- function(parameters, nSites, nLoci, nDip, mutrate, extra = FALSE) {

  # this function is intended to be used with a four-population model
  nPops <- 4

  # Read the vector with the parameters and assign each parameter to the correct command name
  Ne <- parameters[1]
  # set the relative size of each population - for the extant populations
  N1 <- parameters[2]; N2 <- parameters[3]; N3 <- parameters[4]; N4 <- parameters[5]
  # and the ancient populations
  NA1 <- parameters[6]; NA2 <- parameters[7]

  # get the proportion of loci with migration - no barriers against migration between the different ecotypes
  pM <- parameters[19]
  # get the proportion of loci without migration - from the crab to the wave ecotype at the same location
  pCW <- parameters[20]
  # get the proportion of loci without migration - from the wave to the crab ecotype at the same location
  pWC <- parameters[21]
  # and the proportion of loci without any migration - total barrier against migration between the different ecotypes
  pNO <- parameters[22]

  # get the migration rates - between different ecotypes at the same site - this is the m value on the M = 4N0m formula
  # from the crab ecotype to the wave ecotype - at the first site
  mCW1 <- parameters[12]
  # from the crab ecotype to the wave ecotype - at the second site
  mCW2 <- parameters[13]
  # from the wave ecotype to the crab ecotype - at the fist site
  mWC1 <- parameters[14]
  # from the wave ecotype to the crab ecotype - at the second site
  mWC2 <- parameters[15]

  # between crab populations inhabiting different locations
  mCC <- parameters[16]
  # between wave populations inhabiting different locations
  mWW <- parameters[17]
  # and between the two ancestral populations
  mAA <- parameters[18]

  # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
  # and REMEMBER that M = 4N0m
  # set the migration rates to the scale of Nref - between different ecotypes at the same site
  # from the crab ecotype to the wave ecotype - at the first site
  mig_CW1 <- 4*Ne*mCW1 # -m 2 1 mig_CW
  # from the crab ecotype to the wave ecotype - at the second site
  mig_CW2 <- 4*Ne*mCW2 # -m 4 3 mig_CW
  # from the wave ecotype to the crab ecotype - at the first site
  mig_WC1 <- 4*Ne*mWC1 # -m 1 2 mig_WC
  # from the wave ecotype to the crab ecotype - at the second site
  mig_WC2 <- 4*Ne*mWC2 # -m 3 4 mig_WC

  # between crab populations at different locations
  mig_CC <- 4*Ne*mCC
  # between wave populations at different locations
  mig_WW <- 4*Ne*mWW
  # and between the two ancestral populations
  mig_AA <- 4*Ne*mAA

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
  # a variable to inform when does migration start between ancestral populations
  tmAA <- Rsplit + 0.0001
  # if Rsplit is zero, we can not subtract something from it
  if(Rsplit != 0) {
    # when Rsplit is not zero, set the end of the migration before the split
    tmRS <- Rsplit - 0.0001
  } else {
    # when Rsplit is zero, set the end of the migration at the split
    tmRS <- Rsplit
  }

  # use a multinomial distribution to get the number of loci simulated under each category
  lociTotal <- as.vector(stats::rmultinom(n = 1, size = nLoci, c(pM, pCW, pWC, pNO)))

  # if extra is TRUE, then more loci than required per category will be simulated
  if(extra == TRUE) {
    # save the required number of simulated loci per category
    targetLoci <- lociTotal
    # simulate more 25 loci per category
    lociTotal <- lociTotal + 25
  }

  # cheat code: pop1 - crab in site 1; pop2 - wave in site 1; pop3 - crab in site 2; pop4 - wave in site 2

  # create command line with no barriers to migration
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and
  # nrep: the number of independent loci that will be produced
  with.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[1], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                    # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                    "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                    # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
                    # set the migration rate between different ecotypes inhabiting the same location
                    "-m 2 1", mig_CW1, "-m 4 3", mig_CW2, "-m 1 2", mig_WC1, "-m 3 4", mig_WC2,
                    # set the migration rate between the same ecotypes inhabiting different locations
                    "-m 3 1", mig_CC, "-m 1 3", mig_CC, "-m 4 2", mig_WW, "-m 2 4", mig_WW,
                    # set the migration, between all populations, to zero - immediately before the recent split
                    "-eM", tmRS, "0",
                    # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                    "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                    # set the size of the ancient populations
                    # -en <t> <i> <n>: Set the size of population i to n*N0 at time t.
                    "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                    # set the migration rate between the ancestral populations
                    "-em", tmAA, "2 3", mig_AA, "-em", tmAA, "3 2", mig_AA,
                    # set the migration, between all populations, to zero - immediately at the ancient split
                    "-eM", Asplit, "0",
                    # add a split event - this event creates the two ancestral populations
                    "-ej", Asplit, "2 3",
                    # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                    # -eN <t> <n>: set the size of all populations to n*N0 at time t.
                    "-eN", Asplit, 1)

  # create command line with no migration from the crab to the wave ecotype within the same location
  no.mCW <- paste(paste(nDip*2, collapse = " "), lociTotal[2], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                  # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0
                  "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                  # set the migration rate from wave to crab (looking forward in time) inhabiting the same location
                  "-m 2 1 0 -m 4 3 0 -m 1 2", mig_WC1, "-m 3 4", mig_WC2,
                  # set the migration rate between the same ecotypes inhabiting different locations
                  "-m 3 1", mig_CC, "-m 1 3", mig_CC, "-m 4 2", mig_WW, "-m 2 4", mig_WW,
                  # set the migration, between all populations, to zero - immediately before the recent split
                  "-eM", tmRS, "0",
                  # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                  "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                  # set the size of the ancient populations
                  "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                  # set the migration rate between the ancestral populations
                  "-em", tmAA, "2 3", mig_AA, "-em", tmAA, "3 2", mig_AA,
                  # set the migration, between all populations, to zero - immediately at the ancient split
                  "-eM", Asplit, "0",
                  # add a split event - this event creates the two ancestral populations
                  "-ej", Asplit, "2 3",
                  # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                  "-eN", Asplit, 1)

  # create command line with no migration from the wave to the crab ecotype within the same location
  no.mWC <- paste(paste(nDip*2, collapse = " "), lociTotal[3], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                  # set the size of the present day populations
                  "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                  # set the migration rate from crab to wave (looking forward in time) inhabiting the same location
                  "-m 2 1", mig_CW1, "-m 4 3", mig_CW2, "-m 1 2 0 -m 3 4 0",
                  # set the migration rate between the same ecotypes inhabiting different locations
                  "-m 3 1", mig_CC, "-m 1 3", mig_CC, "-m 4 2", mig_WW, "-m 2 4", mig_WW,
                  # set the migration, between all populations, to zero - immediately before the recent split
                  "-eM", tmRS, "0",
                  # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                  "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                  # set the size of the ancient populations
                  "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                  # set the migration rate between the ancestral populations
                  "-em", tmAA, "2 3", mig_AA, "-em", tmAA, "3 2", mig_AA,
                  # set the migration, between all populations, to zero - immediately at the ancient split
                  "-eM", Asplit, "0",
                  # add a split event - this event creates the two ancestral populations
                  "-ej", Asplit, "2 3",
                  # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                  "-eN", Asplit, 1)

  # create a command line for the loci without any migration (between the different ecotypes)
  without.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[4], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "), "0",
                       # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                       "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                       # set the migration rate between the same ecotypes inhabiting different locations
                       "-m 3 1", mig_CC, "-m 1 3", mig_CC, "-m 4 2", mig_WW, "-m 2 4", mig_WW,
                       # set the migration, between all populations, to zero - immediately before the recent split
                       "-eM", tmRS, "0",
                       # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                       "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                       # set the size of the ancient populations
                       "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                       # set the migration rate between the ancestral populations
                       "-em", tmAA, "2 3", mig_AA, "-em", tmAA, "3 2", mig_AA,
                       # set the migration, between all populations, to zero - immediately at the ancient split
                       "-eM", Asplit, "0",
                       # add a split event - this event creates the two ancestral populations
                       "-ej", Asplit, "2 3",
                       # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                       "-eN", Asplit, 1)

  # combine the four different types of commands
  cmdParallel <- c(with.mig, no.mCW, no.mWC, without.mig)

  # if extra is equal to TRUE, then we simulated more loci than required
  if(extra == TRUE)
    # include the required number of loci per category in the output
    cmdParallel <- list(commands = cmdParallel, targetLoci = targetLoci)

  # output the command line for the parallel origin model
  cmdParallel
}

#' Organize scrm output
#'
#' This function is utilized to sort out the scrm output. The order of the
#' populations changes accordingly to the model used (i.e. single or parallel
#' origin). Running this function will re-organize the output produced by scrm,
#' so that the populations are in the same order in both models.
#'
#' @param seg_sites a matrix of segregating sites as produced by scrm. Each
#'   column of the matrix is a different site and each row is a different
#'   haplotype.
#' @param nHap an integer representing the total number of haplotypes simulated.
#' @param nPops an integer, representing the total number of populations of the
#'   simulated model.
#'
#' @return a matrix of segregating sites, similar to `seg_sites` but with the
#'   populations organized so that the order is always the same, regardless of
#'   the model used.
#'
#' @keywords internal
#'
#' @export
organizeSCRM <- function(seg_sites, nHap, nPops) {

  # get the number of haplotypes simulated by population
  haPop <- nHap/nPops
  # create a vector with the index representing the beginning of each population
  beginPop <- seq(from = 1, to = nHap, by = haPop)

  # remove the name (position) of each site - this is something that scrm creates
  seg_sites <- unname(seg_sites)

  # in the single model, we need to switch the order of the second and third population
  # get the haplotypes corresponding to each population
  pop2 <- seg_sites[(beginPop[2]):(beginPop[3]-1), ]
  pop3 <- seg_sites[(beginPop[3]):(beginPop[4]-1), ]

  # re-organize the matrix of haplotypes with the populations in the correct order
  seg_sites[(beginPop[2]):(beginPop[3]-1), ] <- pop3
  seg_sites[(beginPop[3]):(beginPop[4]-1), ] <- pop2

  # output the matrix of haplotypes with the populations in the correct order
  seg_sites
}


#' Run scrm and obtain genotypes
#'
#' This function will run the scrm package, according to the command line
#' supplied as input. It will also combine haplotypes into genotypes and
#' re-organize the output if the simulations were performed under a single
#' origin scenario. This is to ensure that the output of the four-population
#' models will always follow the same order: the two divergent ecotypes in the
#' first location, followed by the two divergent ecotypes in the second
#' location.
#'
#' @param commands A character string containing the commands for the scrm
#'   package. This string can be created using the `cmd2pops`, the `cmdSingle`
#'   or the `cmdParallel` functions.
#' @param nDip An integer representing the total number of diploid individuals
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this.
#' @param nPops An integer that informs of how many populations exist on the
#'   model you are trying to run.
#' @param model Either "2pops", "Single" or "Parallel" indicating which model
#'   should be simulated.
#'
#' @return a list with the simulated genotypes. Each entry is a different locus
#'   and, for each locus, different rows represent different individuals and
#'   each column is a different site.
#'
#' @examples
#'
#' @export
runSCRM <- function(commands, nDip, nPops, model) {

  # check if the input is correct - the selected model should be one of the following
  if(model %in% c("2pops", "Single", "Parallel") == FALSE)
    stop(paste("The selected model should be either 2pops, Single or Parallel. Please check"))

  # binding the variable locally to the function
  temp1 <- NULL

  # run the scrm package for each set of commands - with and without migration
  simulation <- lapply(commands, FUN = function(x) scrm::scrm(x))

  # extract the information from each simulation and store it on a temporary matrix
  for (i in 1:length(simulation)) {

    # create the temporary matrix for each simulation
    assign(paste("temp", i, sep = ""), simulation[[i]][["seg_sites"]])
  }

  if(length(simulation) != 1) {

    # combine all simulations into a matrix of haplotypes
    haplotypes <- append(temp1, unlist(mget(paste0("temp", 2:length(simulation))), recursive = FALSE, use.names = FALSE))

  } else {

    # if only set of simulations was performed, only one set of haplotypes exist
    haplotypes <- temp1
  }

  # get the total number of haplotypes
  nHap <- nDip*2

  # apply a correction for the situations where scrm does not produce a single polymorphic site
  # first check the dimensions of each list entry
  size <- matrix(unlist(lapply(haplotypes, dim)), ncol = 2, byrow = TRUE)

  # if one entry has no columns, i.e. no sites, then add columns containing only zeros to that entry
  if(any(size[, 2] == 0)) {

    # add two columns containing zeros to that locus
    haplotypes <- haplo.fix(haplotypes = haplotypes, nHap = nHap)
  }

  # re-organize output for the single model
  if (model == "Single") {

    # change the order of the populations in the single origin model
    # so that the order is always: ecotype C and W in the first location and ecotype C and W in the second location
    haplotypes <- lapply(haplotypes, function(segSites) organizeSCRM(segSites, nHap, nPops))
  }

  # convert the haplotypes to genotypes
  genotypes <- GetGenotypes(haplotypes, nDip = nDip)

  # output the genotypes
  genotypes
}


#' Force the simulations to contain the required number of loci
#'
#' This function attempts to force the required number of loci after the
#' filtering steps are performed.
#'
#' This is done by simulating extra loci for each of the different types of
#' simulations performed. The possible types of simulations include loci without
#' barriers against migration between divergent ecotypes, loci without migration
#' from the C towards the W ecotype, loci without migration from the W towards
#' the C ecotypes and loci where no migration occurs between divergent ecotypes.
#' Using this function, more loci than required are simulated for each of those
#' types of simulations.
#'
#' Then, a coverage-based filter is applied to the data, followed by a filter
#' based on a required number of minor-allele reads per site. Those filters
#' remove some loci from the data. The extra simulated loci should allow us to
#' keep the required number of loci per type of simulation even after filtering.
#'
#'
#' @param model a character, either 2pops", "Single" or "Parallel" indicating
#'   which model should be simulated.
#' @param parameters a vector of parameters used to create the command line for
#'   the scrm package. Each entry of the vector is a different parameter. Note
#'   that each vector entry should be named with the name of the corresponding
#'   parameter. The output of the `CreateParameters` function is the intended
#'   input.
#' @param nSites is an integer that specifies how many base pairs should scrm
#'   simulate, i.e. how many sites per locus to simulate.
#' @param nLoci an integer that represents how many independent loci should be
#'   simulated.
#' @param nDip an integer representing the total number of diploid individuals
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this. Also note that this is the total
#'   number of diploid individuals and this function will distribute the
#'   individuals equally by the simulated populations.
#' @param mutrate an integer representing the mutation rate assumed for the
#'   simulations.
#' @param mean an integer or a vector defining the mean value of the negative
#'   binomial distribution from which different number of reads are drawn. It
#'   represents the mean coverage across all sites. If a vector is supplied, the
#'   function assumes that each entry of the vector is the mean for a different
#'   population.
#' @param variance an integer or a vector defining the variance of the negative
#'   binomial distribution from which different number of reads are drawn. It
#'   represents the variance of the total coverage across all sites. If a vector
#'   is supplied, the function assumes that each entry of the vector is the
#'   variance for a different population.
#' @param minimum an integer representing the minimum coverage allowed. Sites
#'   where any population has a depth of coverage below this threshold are
#'   removed from the data.
#' @param maximum an integer representing the maximum coverage allowed. Sites
#'   where any population has a depth of coverage above this threshold are
#'   removed from the data.
#' @param size a list with one entry per population. Each entry should be a
#'   vector containing the size (in number of diploid individuals) of each pool.
#'   Thus, if a population was sequenced using a single pool, the vector should
#'   contain only one entry. If a population was sequenced using two pools, each
#'   with 10 individuals, this vector should contain two entries and both will
#'   be 10.
#' @param min.minor is an integer representing the minimum allowed number of
#'   minor-allele reads. Sites that, across all populations, have less
#'   minor-allele reads than this threshold will be removed from the data.
#'
#' @return a list with two names entries
#'
#'   \item{pool}{a list with three different entries: major, minor and total.
#'   This list is obtained by running the \code{\link{forcePool}} function.}
#'
#'   \item{nPoly}{a numeric value indicating the mean number of polymorphic
#'   sites across all simulated locus.}
#'
#' @examples
#'
#'
#' @export
forceLocus <- function(model, parameters, nSites, nLoci, nDip, mutrate, mean, variance, minimum, maximum, size, min.minor) {

  # create the command line to run the scrm package
  # the command line varies according to the selected model
  if(model == "2pops") {

    # create the command line for a star shaped model
    commands <- cmd2pops(parameters, nSites, nLoci, nDip, mutrate, extra = TRUE)
    # set the number of populations
    nPops <- 2

  } else if (model == "Parallel") {

    # create the command line for the parallel origin model
    commands <- cmdParallel(parameters, nSites, nLoci, nDip, mutrate, extra = TRUE)
    # set the number of populations
    nPops <- 4

  } else if (model == "Single") {

    # create the command line for the single origin model
    commands <- cmdSingle(parameters, nSites, nLoci, nDip, mutrate, extra = TRUE)
    # set the number of populations
    nPops <- 4

  } else {

    # if a correct model is not supplied as input for the function - stop and warn
    stop(paste("model should be 2pops, Parallel or Single. Please check!"))
  }

  # get the required number of loci per category - this is the number of loci we want in the end
  target <- commands[["targetLoci"]]
  # get the command line for scrm
  commands <- commands[["commands"]]

  # the length of the commands object indicates how many different categories of simulations we are performing
  # the possible categories are: loci with no barriers to migration, loci with no migration from one ecotype to the other
  # and loci without any migration between the different ecotypes
  nSims <- length(commands)

  # run the scrm package and obtain the genotypes
  genotypes <- lapply(commands, function(sim) runSCRM(commands = sim, nDip, nPops, model))

  # get the mean number of polymorphic sites
  nPoly <- mean(unlist(sapply(genotypes, function(sim) sapply(sim, ncol))))

  # simulate total number of reads per site
  reads <- lapply(genotypes, function(sim) simulateCoverage(mean, variance, genotypes = sim))

  # remove sites with a depth of coverage above or below the defined threshold
  reads <- lapply(1:nSims, function(sim)
    remove_by_reads(nLoci = length(reads[[sim]]), reads[[sim]], minimum, maximum, genotypes = genotypes[[sim]]))

  # get the genotypes - without sites simulated with a coverage below or above the threshold
  genotypes <- lapply(reads, FUN = function(sim) lapply(sim, "[[", 2))
  # get the reads - without sites simulated with a coverage below or above the threshold
  reads <- lapply(reads, FUN = function(sim) lapply(sim, "[[", 1))

  # check the dimensions of the matrices with the genotypes
  # it is possible that some loci do not have any polymorphic site after the coverage filter
  dimensions <- lapply(genotypes, function(sim) sapply(sim, ncol))
  # keep only those loci where we have at least one polymorphic site
  tokeep <- lapply(dimensions, function(sim) sim != 0)
  # remove loci without polymorphic sites from the genotypes
  genotypes <- lapply(1:nSims, function(sim) genotypes[[sim]][tokeep[[sim]]])
  # remove loci without polymorphic sites from the matrices with the coverage
  reads <- lapply(1:nSims, function(sim) reads[[sim]][tokeep[[sim]]])

  # simulate individual contribution to the total number of reads
  indContribution <- lapply(1:nSims, FUN = function(sim) lapply(reads[[sim]], function(locus)
    popsReads(list_np = size, coverage = locus, pError = parameters["PoolError"])))

  # simulate the number of reference reads
  reference <- lapply(1:nSims, FUN = function(sim) lapply(1:length(genotypes[[sim]]), function(locus)
    numberReferencePop(genotypes = genotypes[[sim]][[locus]], indContribution = indContribution[[sim]][[locus]],
                       size = size, error = parameters["SeqError"])))

  # simulate pooled sequencing data
  pool <- lapply(1:nSims, function(sim)
    poolPops(nPops, nLoci=length(indContribution[[sim]]), indContribution=indContribution[[sim]], readsReference=reference[[sim]]))

  # define major and minor alleles
  pool <- lapply(1:nSims, function(sim) lapply(1:length(pool[[sim]][["total"]]), function(locus)
    minorPool(reference = pool[[sim]][["reference"]][[locus]], alternative = pool[[sim]][["alternative"]][[locus]],
              coverage = pool[[sim]][["total"]][[locus]], min.minor)))

  # convert the pool list back to the previous format:
  # one entry for major allele, one for minor allele and a final one for total coverage
  pool <- lapply(1:nSims, function(sim)
    list(major = lapply(pool[[sim]], function(locus) locus[["major"]]), minor = lapply(pool[[sim]], function(locus) locus[["minor"]]),
         total = lapply(pool[[sim]], function(locus) locus[["total"]])))

  # remove loci without polymorphic sites and randomly select the required number of loci per category
  pool <- forcePool(nSims, pool = pool, target)

  # output the pooled sequencing data and the number of polymorphic sites prior to filtering
  list(pool = pool, nPoly = nPoly)
}
