# Simulation of phenotype data for common and rare variants

source(file.path(Sys.getenv("R_SCRIPTDIR"),"SimPhenoBase.R"))

doPrepareTest <- function() {
  dd <<- new.env()
  # Global variables:
  #dd$ft <- 0.05 # Threshold for rare variants # TODO
  dd$plink.geno <- "W:\\data\\work\\iz12\\pipeline\\simulation\\sim1"
  dd$plink.freq <- "W:\\data\\work\\iz12\\pipeline\\simulation\\sim1.simfreq"
  dd$N_cases <- 1000
  dd$N_controls <- 1000
  dd$N <- dd$N_cases + dd$N_controls
  dd$p <- 110
  dd$c <- 0.4
  dd$alpha0 <- 0.01
  dd$pheno.bin <- "W:\\data\\work\\iz12\\pipeline\\simulation\\pheno_bin.txt"
  dd$pheno.cont <- "W:\\data\\work\\iz12\\pipeline\\simulation\\pheno_cont.txt"
}


doPrepareRun <- function() {
  dd <<- new.env()
  # Global variables:
  #dd$ft <- 0.05 # Threshold for rare variants # TODO
  
  if(Sys.getenv("PLINK_GENO") != "") {
    dd$plink.geno <- Sys.getenv("PLINK_GENO")
  } else {
    stop("No plink genetics provided.")
  }
  
  if(Sys.getenv("PLINK_FREQ") != "") {
    dd$plink.freq <- Sys.getenv("PLINK_FREQ")
  } else {
    stop("No plink frequencies provided.")
  }
  
  if(Sys.getenv("PHENO_BIN") != "") {
    dd$pheno.bin <- Sys.getenv("PHENO_BIN")
  } else {
    stop("No output filename for pheno.bin provided.")
  }
  
    
  if(Sys.getenv("PHENO_CONT") != "") {
    dd$pheno.cont <- Sys.getenv("PHENO_CONT")
  } else {
    stop("No output filename for pheno.cont provided.")
  }
  
  if(Sys.getenv("PHENO_COVARIATES") != "") {
    dd$pheno.covariates <- Sys.getenv("PHENO_COVARIATES")
  } else {
    stop("No output files for pheno.covariates provided.")
  }
  
  if(Sys.getenv("NUM_CASES") != "") {
    dd$N_cases <- as.numeric(Sys.getenv("NUM_CASES"))
  } else {
    dd$N_cases <- 1000
  }
  
  if(Sys.getenv("NUM_CONTROLS") != "") {
    dd$N_controls <- as.numeric(Sys.getenv("NUM_CONTROLS"))
  } else {
    dd$N_controls <- 1000
  }
  
  dd$N <- dd$N_cases + dd$N_controls
  
  if(Sys.getenv("N_SNP_CASE") != "") {
    dd$p.case <- as.numeric(Sys.getenv("N_SNP_CASE"))
  } else {
    dd$p.case <- 100
  }
  
  if(Sys.getenv("N_SNP_CONTROL") != "") {
    dd$p.control <- as.numeric(Sys.getenv("N_SNP_CONTROL"))
  } else {
    dd$p.control <- 10
  }
  
  dd$p <- dd$p.case + dd$p.control
  
  if(Sys.getenv("C_BETTA") != "") {
    dd$c <- as.numeric(Sys.getenv("C_BETTA"))
  } else {
    dd$c <- 0.4
  }
  
  if(Sys.getenv("ALPHA_0") != "") {
    dd$alpha0 <- as.numeric(Sys.getenv("ALPHA_0"))
  } else {
    dd$alpha0 <- 0.01
  }
  
}

doRun <- function() { #{{{
  #doPrepareTest()
  doPrepareRun()
  runSimulation(dd)
  saveAll(dd)
}#}}}

