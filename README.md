cophesim 1.1.2
Simulation script to simulate phenotypes for common and rare variants.
© Ilya Y. Zhbannikov, 2016-06-29

##Requirements:
1. Python 2.7
2. Plink http://pngu.mgh.harvard.edu/~purcell/plink/simulate.shtml
3. plinkio https://pypi.python.org/pypi/plinkio

## Continuous trait

Continuous trait regression: y = b0 + X1 + X2 + ... + Xi + betta1*G1 + betta2*G2 + ... + betta_p*Gp + epsilon

Where:
b0 : baseline mean

X1, X2, ..., Xi : covariates (normally-distributed, with mean m and sd = sqrt(m/10))

G1, G2, etc.: values of genetic markers (0, 1, 2) 

betta1, betta2, etc. : regression coefficients, betta_i=cbetta*|log10(MAF)| (MAF - minor allele frequency)

epsilon: normally-distributed residual

## Dichotomous trait

Directly taken from Plink for now.

## Output:

<ouput prefix>_pheno_cont.txt, <ouput prefix>_pheno_bin.txt, <ouput prefix>_id.txt <ouput prefix>_covariates.txt, and plink files: .bed, .bim, .fam, .simfreq

## Usage
```
python cophesim.py -out <ouput prefix> [other parameters]
```

###Other parameters:

ncases : Number of cases

ncontrols : Number of controls

snpn : Number of null SNPs

snpd : Number of disease SNPs

b0 : baseline mean (for continuous phenotype)

cbetta : parameter to calculate betta_i

cov : mean values of covariates, separated by commas (no spaces!)

## Examples

```
python cophesim.py -out ~/Projects/cophesim/test.output -b0 1200

python cophesim.py -out ~/Projects/cophesim/test.output -b0 1200 -cbetta 0.5

python cophesim.py -out ~/Projects/cophesim/test.output -b0 1200 -cbetta 0.5 -cov 56.7,1.23
```
