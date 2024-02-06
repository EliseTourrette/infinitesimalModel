# 28/10/2020
# updated 27-10-2020 to ensure that the FFTs are done only on lattices that are
# powers of 2
# To parallelize the code, could separate out even and odd for the part
# of the code that exploits the gamma and delta coefficients, but one would also
# require parallelisation of the FFTs.
#
# Code based on following not the variances but (1 - variances)
# The iteration equations are the same mathematically, only the boundary 
# estimates and the resetting at the origin change.
# Thus interpret in all this code "sigma2" as being in fact "1-variance"
# This has advantages for numerical accuracy in the far away regions and also
# allows one to easily get the first (1:n^2) and second (1/n^4) analytic 
# corrections to the asymptotic (null) value.

# Note: sum_{n=0}^{infinity} 1/(2n+1)**2 = pi**2/8 = 1.233701
# If we sum over n also negative, the result is double, that is pi**2/4
# The iteration of the variances uses the kernal 
#   K(n,m) = 2/pi^2 1/(distance(n,m))^2  (distance odd)
# because only half of the matter diffuses out, the other half stays put.

# In this code, we use a first lattice (for the variances) ranging from -N 
# to +N (2N+1 sites) that is symmetric and contains the
# (complex Gaussian) variances of the k'th Fourier mode of the genetic value.
# In the mathematical iteration, we need to calculate for each n in [-N,+N]
#   sum_m K(n,m) var(m)
# where m ranges from -infinity to +infinity. Numerically, this sum is 
# broken down into two parts: the one where the site "m" is on the lattice [-N,+N]
# and one where it is beyond. For that second case, we replace var(m) by its 
# approximation (terms gamma/m^2 and delta/m^4) and the sum is then just 
# gamma times a look-up plus delta times another look-up.
# The sum where m belongs to [-N,+N] is done via FFT and requires that one
# extend the lattice so that K(n,m) does not pick up spurious terms. This lattice
# extension is of size 4N+4 (because 2(2N+1) won't work) and it is that lattice 
# size that must be a power of 2 to have a fast FFT calculation.
# 

# R CMD BATCH --no-save --no-restore '--args N=511 log2Generations=17' recurrence_variances.r out.out &

args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  N <- 524287
  log2Generations <- 21 
  alpha <- 1.6
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


build_reference_contribution_1_over_m2_and_m4 <- function(N,CUTOFF) {
  # reference_contribution_from_outside_grid is up to a factor the part of the 
  # infinite sum coming from outside the grid.
  # For any n on the grid (n in [-N, N]) it is defined as 
  # Sum_m K(n,m)/ (m^2 or m^4) for m outside of the grid and where
  # K(n,m) is the kernal = 2/pi^2 1/(distance(n,m))^2 if this distance is odd, 0 otherwise

  # Here it seems that we could store the result just for the sites to the right of 
  # the origin to save a factor 2 is memory space.
  # Code verifications: m2 same as in older code, m4 in agreement with direct calculation when n=0

  reference_contribution_from_outside_grid_1_over_m2 <- rep(NA, 2*N+1)
  reference_contribution_from_outside_grid_1_over_m4 <- rep(NA, 2*N+1)
  for(n in 0:N) {
    vectplus <- rep(1,CUTOFF)
    vectplus <- vectplus/ ((N-n+1):(N+CUTOFF-n))^2    # this is the K(n,m) factor
    vectplus2 <- vectplus/ ((N+1):(N+CUTOFF))^2       # this is the 1/m^2 factor  
    vectplus4 <- vectplus/ ((N+1):(N+CUTOFF))^4       # this is the 1/m^4 factor  
    vectmoins <- rep(1,CUTOFF)
    vectmoins <- vectmoins/ ((N+n+1):(N+CUTOFF+n))^2
    vectmoins2 <- vectmoins/ ((N+1):(N+CUTOFF))^2
    vectmoins4 <- vectmoins/ ((N+1):(N+CUTOFF))^4

    if(n %% 2 == 0) {
      if(N %% 2 == 0) {
        contribution2 <- sum(vectplus2[seq(from = 1, to = CUTOFF, by = 2)]) + sum(vectmoins2[seq(from = 1, to = CUTOFF, by = 2)]) 
        contribution4 <- sum(vectplus4[seq(from = 1, to = CUTOFF, by = 2)]) + sum(vectmoins4[seq(from = 1, to = CUTOFF, by = 2)]) 
      }
      if(N %% 2 != 0) {
        contribution2 <- sum(vectplus2[seq(from = 2, to = CUTOFF, by = 2)]) + sum(vectmoins2[seq(from = 2, to = CUTOFF, by = 2)]) 
        contribution4 <- sum(vectplus4[seq(from = 2, to = CUTOFF, by = 2)]) + sum(vectmoins4[seq(from = 2, to = CUTOFF, by = 2)]) 
      }
    }
    if(n %% 2 != 0) {
      if(N %% 2 == 0) {
        contribution2 <- sum(vectplus2[seq(from = 2, to = CUTOFF, by = 2)]) + sum(vectmoins2[seq(from = 2, to = CUTOFF, by = 2)]) 
        contribution4 <- sum(vectplus4[seq(from = 2, to = CUTOFF, by = 2)]) + sum(vectmoins4[seq(from = 2, to = CUTOFF, by = 2)]) 
      }
      if(N %% 2 != 0) {
        contribution2 <- sum(vectplus2[seq(from = 1, to = CUTOFF, by = 2)]) + sum(vectmoins2[seq(from = 1, to = CUTOFF, by = 2)]) 
        contribution4 <- sum(vectplus4[seq(from = 1, to = CUTOFF, by = 2)]) + sum(vectmoins4[seq(from = 1, to = CUTOFF, by = 2)]) 
      }
    }
    reference_contribution_from_outside_grid_1_over_m2[N+1+n] <- contribution2
    reference_contribution_from_outside_grid_1_over_m4[N+1+n] <- contribution4
  }
  reference_contribution_from_outside_grid_1_over_m2[1:N] <- reference_contribution_from_outside_grid_1_over_m2[(2*N+1):(N+2)]
  reference_contribution_from_outside_grid_1_over_m4[1:N] <- reference_contribution_from_outside_grid_1_over_m4[(2*N+1):(N+2)]
  return(list((2/pi^2)*reference_contribution_from_outside_grid_1_over_m2, (2/pi^2)*reference_contribution_from_outside_grid_1_over_m4))
}  

get_gammas_deltas <- function(N,sigma2) {
  # Estimation of the parameters gamma and delta for the profile passed.
  # In the tails, we assume the profile behaves as gamma/m^2 + delta/m^4, separately 
  # for even and odd sites. The two gammas and two deltas are estimated via
  # the sites at N and 7N/8 to the left of the origin, based on the system
  #  ( 1/d1^2    1/d1^4 ) ( gamma )    =   ( sigma2(site1)  ) 
  #  ( 1/d2^2    1/d2^4 ) ( delta )    =   ( sigma2(site2)  ) 
  # so that sigma2 at each site is asymptotically the sum of a part in 1/m^2 and a part in 1/m^4
  # case of even sublattice
  if(N %% 2 == 0) {
    siteN <- 1
    site7Nover8 <- 2*round(N/16)+1
  }
  if(N %% 2 != 0) {
    siteN <- 2
    site7Nover8 <- 2*round(N/16)
  } 
  d1 <- N+1-siteN
  d2 <- N+1-site7Nover8
  determinant <- (1/d2^2 - 1/d1^2) / (d1*d1*d2*d2)
  gamma_even <- (sigma2[siteN]/d2^4 - sigma2[site7Nover8]/d1^4)/determinant
  delta_even <- (-sigma2[siteN]/d2^2 + sigma2[site7Nover8]/d1^2)/determinant
  # case of odd sublattice
  siteN <- siteN+1
  site7Nover8 <- site7Nover8 + 1
  d1 <- N+1-siteN
  d2 <- N+1-site7Nover8
  determinant <- (1/d2^2 - 1/d1^2) / (d1*d1*d2*d2)
  gamma_odd <- (sigma2[siteN]/d2^4 - sigma2[site7Nover8]/d1^4)/determinant
  delta_odd <- (-sigma2[siteN]/d2^2 + sigma2[site7Nover8]/d1^2)/determinant

  return(list(gamma_even,gamma_odd,delta_even,delta_odd))
}

selectionRound <- function(N, mu, alpha, sigma2, reference_contribution_from_outside_grid_1_over_m2, reference_contribution_from_outside_grid_1_over_m4, TdiffusionFT, extendedSigma2) {  
  # reference_contribution_from_outside_grid is up to a factor the part of the infinite sum coming from outside the grid
  # For any n on the grid (n in [-N, N]) it is defined as Sum_m K(n,m)/ (m^2 or m^4) for m outside of the grid and
  # where K(n,m) is the kernal = 2/pi^2 1/(distance(n,m))^2 if this distance is odd, 0 otherwise 
  
  newSigma2 <- rep(NA, 2*N+1)               # Ultimately it may be good to not have that factor 2 of redundancy
                                            # but since in practice it is not memory that is limiting probably
                                            # there is no point in reducing this redundancy
  
  # selection
  mu <- mu + alpha*sqrt(sigma2[N+1])        # Select for individuals that are alpha (current) standard deviations better than the current mean

#sigma2 <- rep(0, 2*N+1) 

  sigma2[N+1] <- 1                          # Using that selection rather than truncation selection keeps all the distributions Gaussian
                                            # Note this Dirac delta function selection can be interpreted as a Gaussian selection with infinitesimal width
  

  # In this function, apply one iteration of the recurrence on the quantities (1-variances here)      
                    
  # First deal with the sites present on the lattice
  # To do the convolution within the sites kept explicitly, we resort to FFT
  # The current set of sites is embedded in a lattice of double the size plus 2 extra sites
  # so that we have 2(2N+1)+2 sites in total, 
  # that is 2^p = 4*(N+1) with periodic boundary conditions
  # On these additional sites we set the variances to 0 so that the
  # convolution on the whole extended lattice gives exactly the convolution
  # of the internal sites
  
  extendedSigma2[(N+2):(N+1+2*N+1)] <- sigma2

    
  extendedSigma2FT <- fft(extendedSigma2, inverse = FALSE)
  
  vectProduct <- extendedSigma2FT * TdiffusionFT
  
  extendedSigma2 <- Re(fft(vectProduct, inverse = TRUE)/(2*(2*N+2)))
  newSigma2 <- 0.5*sigma2 + extendedSigma2[(N+2):(N+2*N+2)]



  # Moving the next part out of the previous loop hardly improves the speed, clearly
  # cpu time is dominated but the sum above (and numbers are exactly the same)

  # Second deal with the contributions from the sites not present on the lattice
  # Need to get the fits to tails of the form gamma/m^2 + delta/m^4, separately 
  # for even and odd sites.
  gammas_deltas <- get_gammas_deltas(N,sigma2)
  gamma_even <- gammas_deltas[[1]]
  gamma_odd  <- gammas_deltas[[2]]
  delta_even <- gammas_deltas[[3]]
  delta_odd  <- gammas_deltas[[4]]

  RightSeqImpair <- N+1+seq(from = 1, to = N, by = 2)
  RightSeqPair <- N+1+seq(from = 0, to = N, by = 2)
  newSigma2[RightSeqImpair] <- newSigma2[RightSeqImpair] + gamma_even*reference_contribution_from_outside_grid_1_over_m2[RightSeqImpair] +
delta_even*reference_contribution_from_outside_grid_1_over_m4[RightSeqImpair]
  newSigma2[RightSeqPair] <- newSigma2[RightSeqPair] + 
gamma_odd*reference_contribution_from_outside_grid_1_over_m2[RightSeqPair] +
delta_odd*reference_contribution_from_outside_grid_1_over_m4[RightSeqPair]
  # Note: for constructing g=1 profile, gamma=delta=0
  # for constructing g=2 profile, gamma_odd = 0.2026423, others = 0
  # for constructing g=3 profile, gamma_odd is 2/3 of before, delta_odd=0 (to machine precision)
  # at that generation gamma_even is unchanged because pure Cauchy diffusion leads to its doubling
  # but only half of odd profile participates to diffusion

  newSigma2[1:N] <- newSigma2[(2*N+1):(N+2)]       # Use symmetry to half CPU time

  return(list(mu, newSigma2))
}

########################################################################
########################################################################     
###                  MAIN PART OF THE CODE                          ####
########################################################################
########################################################################

# Probably should take N greater than Generations

##! N <- 511  # Your proposed value of N
# Readjust so that 4(N+1) = 2^pow:
pow <- round(log2(4*(N+1)))
N = 2**(pow-2) - 1  # Readjusted value so that the FFT is efficient
print(N)


reference_contribution_from_outside_grid <- build_reference_contribution_1_over_m2_and_m4(N,100000)
reference_contribution_from_outside_grid_1_over_m2 <- reference_contribution_from_outside_grid[[1]]
reference_contribution_from_outside_grid_1_over_m4 <- reference_contribution_from_outside_grid[[2]]


# These contributions are peaked for "n" at the borders of the grid as anticipated.
# plot(reference_contribution_from_outside_grid_1_over_m2)

sigma2 <- rep(NA, 2*N+1)

# Proposed modified initial conditions to reduce transient behavior
# sigma2[N+1] <- 0.95  # value at origin, in fact this is not used since it
                    # it is reset to 0 at each time step (we use variances here, not 1-variances)
# for all the rest, use a hyperbolic tangent function with scale "scale"
# scale <- 0.01 
# sigma2[(N+2):(2*N+1)] = 1 - (exp(scale*(1:N))-exp(-scale*(1:N)))/(exp(scale*(1:N))+exp(-scale*(1:N)))
# sigma2[1:N] = sigma2[(2*N+1):(N+2)]
# That law does not lead to good results, indeed the curves go down and then back up. Clearly
# the shape at long times is more like exp(-x/L) than a Gaussian: this 1/distance**2 diffusion
# does not lead to a curve with an inflexion point => would have to change this initial shape, e.g.,
# 1 - 1 / (1 *  x**2)

# Can also use instead the usual initial conditions to control the first generations
sigma2 <- rep(0, 2*N+1)   #  0 because here we work with 1 - variance

# alpha <- 1.6  # selection point: we keep only parents that have a value alpha*current_standard_deviation above current mean.
mu <- 0     # initial value of the mean in the distribution of genetic values

##! log2Generations <- 17
Generations <- 2**log2Generations    # total number of iterations to use
listMu <- rep(NA,Generations+1)
g <- Generations

sigmaSquare <- data.frame(matrix(NA, ncol = (1+log2Generations), nrow = N+1))  # only store the right half of grid
names(sigmaSquare) <- 1:(1+log2Generations)      #  column names are just 1 plus the log2 generation number
row.names(sigmaSquare) <- 0:N  #  row names are the index of (n)

gammas_deltas <- data.frame(matrix(NA, ncol = 4, nrow = (1+log2Generations)))  
row.names(gammas_deltas) <- 1:(1+log2Generations)      #  row names are just 1 + the log2 generation number
names(gammas_deltas) <- c("gamma_even","gamma_odd","delta_even","delta_odd")  #  column names 

listSigma0 <- rep(NA, g)

## do the diffusion and its FFT before the for loop -> gain of time
seqImpair <- seq(from = 1, to = 2*N, by = 2)
Tdiffusion <- rep(0, 4*(N+1))
Tdiffusion[1 + seqImpair] <- 1/(seqImpair)^2 
Tdiffusion[2*(2*N+2)+1-seqImpair] <- 1/seqImpair^2
Tdiffusion <- Tdiffusion*2/pi^2
TdiffusionFT <- fft(Tdiffusion, inverse = FALSE)

## avoid memory allocation to extendedSigma2 at every generation
extendedSigma2 <- rep(0,4*(N+1))  # Think of circular lattice having:
# (N+1) sites at 0, then the 2N+1 sites for the sigmas, then again (N+1) sites and finally
# the "image" of site "0" to make the 4(N+1) sites 


time <- proc.time()

for (j in 1:g) {
  generation <- selectionRound(N, mu, alpha, sigma2,reference_contribution_from_outside_grid_1_over_m2,reference_contribution_from_outside_grid_1_over_m4, TdiffusionFT, extendedSigma2)
  mu <- as.vector(generation[[1]])
  sigma2 <- as.vector(generation[[2]])
  listMu[j] <- mu
  listSigma0[j]<-sigma2[N+1]
  # Store and save to disk the profiles if j is a power of 2
  log2gen <- log(j)/log(2)
  if(j == 2**floor(log2gen)) {
    # compute the gammas and deltas for this generation (not the previous one)
    these_gammas_deltas <- get_gammas_deltas(N,sigma2)
    gammas_deltas[1+floor(log2gen),1] <- these_gammas_deltas[[1]]
    gammas_deltas[1+floor(log2gen),2] <- these_gammas_deltas[[2]]
    gammas_deltas[1+floor(log2gen),3] <- these_gammas_deltas[[3]]
    gammas_deltas[1+floor(log2gen),4] <- these_gammas_deltas[[4]]
    sigmaSquare[,1+floor(log2gen)] <- sigma2[(N+1):(2*N+1)]  # we only store the iterations that are powers of 2
    saveRDS(sigmaSquare, file=paste("sigmaSquare2_N",N,"_log2Gen",log2Generations,"_alpha", alpha, ".rds",sep=""))
    write.table(matrix(c(j,mu),nrow=1), file=paste("mu2_N",N,"_log2Gen",log2Generations,"_alpha", alpha, ".txt",sep=""), append = TRUE)
  }
  if(j <= 10000)  write.table(matrix(c(j,mu),nrow = 1), file=paste("mu2_N",N,"_log2Gen",log2Generations,"_alpha", alpha, ".txt",sep=""), append = TRUE)
  if(j <= 10000)  write.table(matrix(c(j,sigma2[N+1]),nrow = 1), file=paste("mu2_N",N,"_log2Gen",log2Generations,"_alpha", alpha, ".txt",sep=""), append = TRUE)
}

cputime <- proc.time() - time



cputime

# Summary of CPU times in seconds (loop over g) for different N and log2gen
#  N \ log2gen   14     15      16      17
#      500       56     116            466
#     1000              405           1559
#     2000             1490
#     3000                      6769
#     4000
#     6000                     25806

# Now using the 4(N+1) being a power of 2
# Summary of times in seconds for different N and log2gen
#  N \ log2gen   14     15      16      17
#      511       3.4    6.9             27.6
#     1023       6.5   13.3             52.5
#     2000             
#     3000                      
#     4000
#     6000                     



