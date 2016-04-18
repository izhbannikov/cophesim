# Simulation script to simulate phenotypes for common and rare variants.
# © Ilya Y. Zhbannikov
# 2016-03-14
#
#--------------------------------------------------------------------------------------------------------#
#                                     Phenotype simulation                                               #
# (for more information see the papers: "Kernel Machine SNP-Set Testing Under Multiple Candidate Kernels"#
#  by Wu et al, 2013)  and "Rare-Variant Association Testing for Sequencing Data with the Sequence Kernel 
# Association Test" by Wu et al, 2011.                                                                   #
#--------------------------------------------------------------------------------------------------------#

# Continuous trait: y = 0.5*X1 + 0.5*X2 + betta1*G1 + betta2*G2 + ... + betta_p*Gp + epsilon
#
# X1, X2 : covariates
# G1, G2, etc.: values of genetic markers (0, 1, 2) 
# betta1, betta2, etc. : regression coefficients (betta0 = 0 for simplicity).betta_i=c*|log10(MAF)|
# epsilon: normally-distributed residual
#
# Dichotomous trait: logit P(y=1) = alpha0 + 0.5*X1 + 0.5*X2 + betta1*G1 + betta2*G2 + ... + betta_p*Gp + epsilon
#
#
#-----------------------------------------Usage-----------------------------------------------------------#
#	1. Change the values of N_cases, N_controls, snp_null, snp_disease, etc. in simulation.sh 
#	2. Run simulation script: ./simulation.sh 
#	3. See output files: pheno_bin.txt, pheno_cont.txt, covariates.txt
#


