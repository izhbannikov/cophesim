#!/bin/bash

SCRIPTDIR=$(readlink -f $(dirname $0))

N_cases=1000
N_controls=1000
snp_null=100
snp_disease=10
cbetta=0.004
alpha0=0.01

echo "Making wgas.sim configuration file..."
rm wgas.sim
echo "$snp_null null 0.00 1.00 1.00 1.00" >> wgas.sim
echo "$snp_disease disease 0.00 1.00 2.00 mult" >> wgas.sim
echo "Done."

# Simulation genotypes with Plink.
# For help visit: http://pngu.mgh.harvard.edu/~purcell/plink/simulate.shtml
plink --simulate wgas.sim --make-bed --out sim1 --simulate-ncases $N_cases --simulate-ncontrols $N_controls

echo "Running SimPheno"

$SCRIPTDIR/./SimPheno $SCRIPTDIR/sim1 $SCRIPTDIR/sim1.simfreq $SCRIPTDIR/pheno_bin.txt $SCRIPTDIR/pheno_cont.txt $SCRIPTDIR/covariates.txt $N_cases $N_controls $snp_null $snp_disease $cbetta $alpha0

echo "Simulation finished."
