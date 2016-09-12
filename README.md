# cophesim

Simulation tool to simulate phenotypes: continuous, dichotomous and survival for common variants from existing genotype data simulated with some other tool.

Version: 1.2.1

© Ilya Y. Zhbannikov, 2016-09-10

##Requirements

1. Python 2.7
2. plinkio https://pypi.python.org/pypi/plinkio

## Input

Files in Plink format: .bed, .bim, .fam.

## Output:

<ouput prefix>_pheno_cont.txt, <ouput prefix>_pheno_bin.txt, <output prefix_surv.txt>, <ouput prefix>_id.txt <ouput prefix>_covariates.txt

## Usage
```
python cophesim.py -i <path to the directory with .bed, .bim, .fam files> -out <ouput prefix> [other parameters]
```

###Other parameters:

-d A flag that indicates dichotomous outcome (on by default)
-c A flag that indicates continuous outcome (off by default)
-s A flag that indicates survival outcome (off by default)


## Examples

```
python cophesim.py -out ~/Projects/cophesim/test.output

python cophesim.py -out ~/Projects/cophesim/test.output -c -s

```
