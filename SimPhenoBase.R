# Simulation script to simulate phenotypes for common and rare variants.
# © Ilya Y. Zhbannikov
# 2016-03-14

library(MultiPhen)


beta_p <- function(c, maf) {
  ans <- c*abs(log10(maf))
  ans
}

sim.covariates <- function(dd) {
  N <- dd$N
  X1 <- runif(N, 0, 1)
  X2 <- rnorm(N, 0, 1)
  cov <- cbind(X1, X2)
  colnames(cov) <- c("X1", "X2")
  dd$covariates <- cov
}

#--------------------------------------------------------------------------------------------------------#
#                                     Phenotype simulation                                               #
# (for more information see the papers: "Kernel Machine SNP-Set Testing Under Multiple Candidate Kernels"#
#  by Wu et al, 2013)  and "Rare-Variant Association Testing for Sequencing Data with the Sequence Kernel 
# Association Test" by Wu et al, 2011.                                                                   #
#--------------------------------------------------------------------------------------------------------#

# Continuous trait: y = 0.5*X1 + 0.5*X2 + betta1*G1 + betta2*G2 + ... + betta_p*Gp + epsilon
#
# X1, X2 : covariates
# betta1, betta2, etc. : regression coefficients (betta0 = 0 for simplicity).betta_i=c*|log10(MAF)|
# epsilon: normally-distributed residual
#
sim.cont.trait <- function(dd) {
  N <- dd$N
  p <- dd$p
  c <- dd$c
  alpha0 <- dd$alpha0
  y.cont <- matrix(nrow=N, ncol=1)
  
  for(i in 1:N) {
    
    X1 <- dd$covariates[i,1]
    X2 <- dd$covariates[i,2]
    
    epsilon <- rnorm(1, 0, 1)
  
    genetics <- sum(unlist(lapply(1:p, function(n) {
                                        beta_p(c, dd$plink.freq.table$V3[n])*dd$plink.data.table[i,n]
                                     })
                                )
                        )

    y.cont[i,1] <- round(0.5*X1 + X2 + genetics + epsilon, 2)
  }
  fam.id <- seq(1,N) # Family id
  dd$y.cont.t <- cbind(fam.id, rownames(dd$plink.data.table), y.cont)
}

# Dichotomous trait: logit P(y=1) = alpha0 + 0.5*X1 + 0.5*X2 + betta1*G1 + betta2*G2 + ... + betta_p*Gp + epsilon
sim.dich.trait <- function(dd) {
  # http://stats.stackexchange.com/questions/46523/how-to-simulate-artificial-data-for-logistic-regression
  N <- dd$N
  p <- dd$p
  c <- dd$c
  alpha0 <- dd$alpha0
  y.dich <- matrix(nrow=N, ncol=1)
  #c <- log(5)/4
  for(i in 1:N) {
    
    X1 <- dd$covariates[i,1]
    X2 <- dd$covariates[i,2]
    
    epsilon <- rnorm(1, 0, 1)
  
    genetics <- sum(unlist(lapply(1:p, function(n) {
                                        beta_p(c, dd$plink.freq.table$V3[n])*dd$plink.data.table[i,n]
                                      })
                        )
                  )
  
    
    z <- alpha0 + 0.5*X1 + 0.5*X2 + genetics + epsilon
    pr <- 1/(1+exp(-z))
    y.dich[i,1] <- rbinom(1,1,pr) 
  }
  fam.id <- seq(1,N) # Family id
  dd$y.dich.t <- cbind(fam.id, rownames(dd$plink.data.table), y.dich)
}

# Saving the data:
saveAll <- function(dd) {
  fam.id <- seq(1,dd$N) # Family id
  write.table(file=dd$pheno.covariates, x=cbind(fam.id, rownames(dd$plink.data.table), dd$covariates), 
              quote = FALSE, row.names=FALSE, col.names = FALSE)
  write.table(file=dd$pheno.cont, x = dd$y.cont.t, 
              quote = FALSE, row.names=FALSE, col.names = FALSE)
  write.table(file=dd$pheno.bin, x = dd$y.dich.t, 
              quote = FALSE, row.names=FALSE, col.names = FALSE)
}

runSimulation <- function(dd){
  # Read plink files:
  dd$plink.data.table <- read.plink(dd$plink.geno) 
  # Reading allele frequencies:
  dd$plink.freq.table <- read.table(dd$plink.freq)
  # Simulate covariates:
  sim.covariates(dd)
  # Simulate continuous trait:
  sim.cont.trait(dd)
  # Simulate dichotomous trait:
  sim.dich.trait(dd)
}

