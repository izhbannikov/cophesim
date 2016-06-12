Simulation script to simulate phenotypes for common and rare variants.
© Ilya Y. Zhbannikov, 2016-06-12

##Requirements:
1. Python 2.7
2. Plink http://pngu.mgh.harvard.edu/~purcell/plink/simulate.shtml
3. plinkio https://pypi.python.org/pypi/plinkio

## Continuous trait

Continuous trait regression: y = b0 + betta1*G1 + betta2*G2 + ... + betta_p*Gp + epsilon

Where:
b0 : baseline mean

G1, G2, etc.: values of genetic markers (0, 1, 2) 

betta1, betta2, etc. : regression coefficients, betta_i=cbetta*|log10(MAF)| (MAF - minor allele frequency)

epsilon: normally-distributed residual

## Dichotomous trait

Directly taken from Plink for now.

## Output:

pheno_cont.txt, pheno_bin.txt, id.txt and plink files: .bed, .bim, .fam, .simfreq

## Usage
```
python cophesim.py -in <path to plink files with file prefix> -out <ouput prefix> [other parameters]
```

###Other parameters:

ncases : Number of cases

ncontrols : Number of controls

snpn : Number of null SNPs

snpd : Number of disease SNPs

b0 : baseline mean (for continuous phenotype)

cbetta : parameter to calculate betta_i

## Examples

python cophesim.py -in ~/Projects/cophesim/simulated.test -out ~/Projects/cophesim/test.output -b0 1200

python cophesim.py -in ~/Projects/cophesim/simulated.test -out ~/Projects/cophesim/test.output -b0 1200 -cbetta 0.5


