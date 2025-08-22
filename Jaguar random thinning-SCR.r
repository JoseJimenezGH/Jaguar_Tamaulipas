#==============================================================================#
#                                                                              #
#    Jaguar density at the northeastern limit of its distribution in México    #
#               RANDOM THINNING-SPATIAL CAPTURE-RECAPTURE                      #
#    Zavdiel A. Manuel-De la Rosa, Leroy Soria-Díaz, Carlos Barriga-Vallejo,   #
#           Gabriela Mendoza-Gutiérrez, Nayeli Martínez-González,              #
#            Claudia C. Astudillo-Sánchez, José Jiménez                        #
#                           18:08 19/04/2025                                   #
#                                                                              #
#==============================================================================#

# Set working directory
setwd('C:/...')

# Load required libraries for analysis
library(scrbook)    # Library for spatial capture-recapture modeling
library(coda)       # Library for MCMC diagnostics
library(mcmcOutput) # Library for MCMC outputs
library(nimble)     # Library for hierarchical modeling
library(terra)      # Spatial data handling

# Load additional functions and scripts
source('Functions/SCR_functions.R')  # Custom SCR-related functions

## SETTING AND PLOTING DATA
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load capture-recapture history data
jaguar.ch <- secr::read.capthist("BData/capt.txt", "BData/traps.txt", detector='count', noccasions=91)

# Summarize capture data and test for population closure
summary(jaguar.ch)
secr::closure.test(jaguar.ch) # Population closure test

# Rearrange capture data for analysis
y3d <- aperm(jaguar.ch, c(1, 3, 2))  # Reorganize capture matrix

# Load trap locations and normalize coordinates
traplocs <- as.matrix(secr::traps(jaguar.ch))
X <- data.matrix(traplocs) / 1000  # Normalize to kilometers
X[,1] <- X[,1] - mean(X[,1])       # Center X-coordinates
X[,2] <- X[,2] - mean(X[,2])       # Center Y-coordinates
rownames(X) <- 1:52                # Rename rows
colnames(X) <- c("X", "Y")         # Set column names

# Sampling information
J <- nrow(X)             # Number of traps
K <- dim(y3d)[3]         # Number of sampling occasions
nind <- dim(y3d)[1]      # Number of detected individuals
detections <- sum(y3d)   # Total number of detections

# Set up data augmentation (M is augmented population size)
M <- 75

# Calculate buffered state space size for irregular trap array
trapShape <- vect("GIS/Traps.shp")
buff_trap <- buffer(trapShape, width = 9150) # 3 * sigma buffer size
buffTrap <- aggregate(buff_trap) # Combine geometries into single polygon
plot(buffTrap)                   # Plot buffered trap array
points(trapShape)                # Add original trap locations to the plot
area <- expanse(buffTrap) / 1e6  # Calculate area in square kilometers

buff <- 9150 / 1000      # Convert buffer to kilometers
xl <- min(X[,1]) - buff  # Calculate x-axis lower bound
xu <- max(X[,1]) + buff  # Calculate x-axis upper bound
yl <- min(X[,2]) - buff  # Calculate y-axis lower bound
yu <- max(X[,2]) + buff  # Calculate y-axis upper bound
xlim <- c(xl, xu)        # Define x-axis limits
ylim <- c(yl, yu)        # Define y-axis limits

# Prepare capture data for augmented individuals
y <- apply(y3d, c(1, 2), sum)  # Sum captures across occasions
yaug <- array(0, c(M, J))      # Augmented capture data
yaug[1:nind, ] <- y            # Fill with observed data

# Visualize capture locations
plot(X, xlim=c(-10, 11), ylim=c(-32, 34), pch="+", cex=0.8, col="red",asp=TRUE)
datn <- apply(y3d, c(2, 3), sum)     # Number of captures at each trap
tot <- apply(datn, 1, sum)           # Total captures by trap
symbols(X, circles=tot/1, inches=F, bg="#00000022", fg=NULL, add=T)
points(X, pch="+", cex=0.8, col="red")# Add trap locations

# Camera trap operation mask
KT <- read.csv("Bdata/traps.csv", sep=",")
KT <- KT[,4:94]  # Columns representing active/inactive status
KT <- data.matrix(KT)  # Convert to matrix
colnames(KT) <- 1:91   # Number of sampling occasions
image(1:K, 1:J, t(KT), yaxt = "n", xlab="Occasion", ylab="", cex.lab=1.25, col=topo.colors(2))
mtext(side = 2, "Camera trap", line = 2.5, cex=1.25)
axis(2, rev(seq(1, J, by=2)))


# Non-ID capture frequencies
jaguar.un <- secr::read.capthist("BData/NID.txt", "BData/traps.txt", detector='count', noccasions=91)
# Load non-ID capture history data with the same number of occasions as ID captures
summary(jaguar.un)
y.un <- aperm(jaguar.un, c(1, 3, 2))  # Rearrange capture data
nnid <- apply(y.un, c(2, 3), sum)     # Non-ID events across traps and occasions
sum(nnid)                             # Total number of non-ID events
nnidd <- apply(nnid, 1, sum)          # Sum of non-ID events by trap


## NIMBLE CODE
##~~~~~~~~~~~~~~
## Define the SCR model
NimModel <- nimbleCode({
  # Parameters
  psi ~ dbeta(1,1)        # Probability of individual inclusion in the population
  lam0 ~ dunif(0,5)       # Baseline detection rate
  sig ~ dunif(0,10)       # Movement parameter
  id.prob ~ dunif(0,1)    # Probability of individual identification

  # Model for augmented individuals
  for(i in 1:M) {
    z[i] ~ dbern(psi)                # Latent inclusion indicator
    s[i,1] ~ dunif(xlim[1], xlim[2]) # x-coordinate of activity center
    s[i,2] ~ dunif(ylim[1], ylim[2]) # y-coordinate of activity center
    
    d2[i,1:J] <- (s[i,1] - X[1:J,1])^2 + (s[i,2] - X[1:J,2])^2
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], 
	                           X = X[1:J,1:2], 
                                   J = J, 
                                   sigma = sig,
                                   lam0 = lam0,
                                   z = z[i])

    # Model for complete capture histories
    y.full[i,1:J] ~ dPoissonVector(lambda = lam[i,1:J] * KT[1:J])
    
    # Observed capture histories
    for(j in 1:J) {
      y.obs[i,j] ~ dbin(id.prob, y.full[i,j])  # Binomial model for ID events
    }

    # Implement zero-trick for irregular trap array
    for(j in 1:J) { # zero-trick by R. Chandler in
      # https://groups.google.com/g/spatialcapturerecapture/c/NzqUovn8jF0/m/Plg2g6O6AgAJ   
      outj[i,j] <- sqrt(d2[i,j]) > buffer  # avoid MCMC sampling beyond buffer distance
    }
    out[i] <- equals(sum(outj[i,1:J]), J)  # Zero-trick condition
    zeros[i] ~ dbern(out[i])               # Zero-trick implementation
  }

  # Derived parameters
  N <- sum(z[1:M])      # Total population size
  D <- 100 * N / area   # Population density (individuals per 100 sq. km)
})

## CONSTANTS
##~~~~~~~~~~~~
KT<-apply(KT,1,sum)
## Constants for the model
str(constants <- list(nnid = nnid,         # Non-ID events
                      J = J,               # Number of traps
                      M = M,               # Data augmentation size
                      KT = KT,             # Sampling occasions
                      area = area          # State space size
))

## DATA
##~~~~~~
## Data for the model
yred <- apply(yaug, c(1, 2), sum)
str(data   <-    list(y.obs = yred,        # ID capture histories
                      buffer = buff,       # Buffer size
                      zeros = rep(0, M),   # Zero-trick vector
                      xlim = xlim,         # State space x-limits
                      ylim = ylim,         # State space y-limits
                      X = X                # Trap coordinates
))

## INITS
##~~~~~~~~
## Initial values for the model
ys <- apply(yaug, c(1, 2), sum)
s.start <- X[sample(1:J, size = M, replace = TRUE), ]
d <- e2dist(s.start[1:M, ], X)   # Euclidean distance between activity centers and traps
lam0s <- runif(1, 0, 0.1)        # Initial baseline detection rate
sigs <- runif(1, 1, 3)           # Initial movement parameter

lam <- lam0s * exp(-d^2 / (2 * sigs^2))

yi <- array(0, c(M, J, K)) # Resighting array
for (j in 1:J) {
  for (k in 1:K) {
    if (nnid[j, k] > 0) {
      probs <- lam[, j] / sum(lam[, j])  # Sampling probabilities
      latent.id <- sample(1:M, nnid[j, k], prob = probs, replace = FALSE)
      yi[latent.id, j, k] <- 1
    }
  }
}

yis <- apply(yi, c(1, 2), sum) + apply(yaug, c(1, 2), sum)
zst <- apply(yis, 1, sum)
zst[zst > 0] <- 1  # Individuals with detections
id.prob.s <- sum(yaug) / (sum(yaug) + sum(nnid))  # Initial ID probability

set.seed(1960)
str(inits   <-   list(z = zst,              # Latent inclusion indicator
                      s = s.start,          # Activity center coordinates
                      lam0 = lam0s,         # Baseline detection rate
                      sig = sigs,           # Movement parameter
                      id.prob = id.prob.s,  # Detection rate for ID events
                      psi = runif(1, 0, 1), # Inclusion probability
                      y.full = yis          # Latent true capture histories
))

## PARAMETERS
##~~~~~~~~~~~~
## Parameters to monitor during MCMC
params <- c('psi', 'lam0', 'sig', 'N', 'D', 'id.prob') # Model parameters

## COMPILATION
##~~~~~~~~~~~~~~
# Compile the model using NIMBLE
start.time <- Sys.time()
Rmodel <- nimbleModel(code = NimModel, constants = constants, data = data, 
                      inits = inits, check = FALSE)
Rmodel$initializeInfo()  # Initialize model information
Rmodel$calculate()       # Calculate initial model state
Cmodel <- compileNimble(Rmodel) # Compile the model

# Configure MCMC settings
conf <- configureMCMC(Rmodel, monitors = params, thin = 10, useConjugacy = TRUE)

# Custom samplers for capture data and activity centers
conf$removeSampler("y.full")
for(j in 1:J) {
  conf$addSampler(target = paste(paste("y.full[1:", M, ", ", j, "]"), sep = ""),
                  type = 'IDSampler', control = list(nnidd = nnidd[j], j = j, M = M), 
                  silent = TRUE)
}

# sampler to jointly update AC's x,y coordinates
conf$removeSampler("s")
ACnodes <- paste0("s[", 1:constants$M, ", 1:2]")
for(node in ACnodes) {
  conf$addSampler(target = node, type = "RW_block", 
                  control = list(adaptScaleOnly = TRUE), silent = TRUE)
}

# Build and compile MCMC object
Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run MCMC sampling
start.time2 <- Sys.time()
outNim <- runMCMC(Cmcmc, niter = 50000, nburnin = 10000, nchains = 3, 
                  inits = inits, setSeed = TRUE, progressBar = TRUE, 
                  samplesAsCodaMCMC = TRUE)
end.time <- Sys.time()
end.time - start.time2  # Time taken for sampling

## OUTPUTS
##~~~~~~~~~~
# Summarize MCMC outputs
summary(mcmcOutput(outNim))
diagPlot(mcmcOutput(outNim))  # Diagnostic plots








