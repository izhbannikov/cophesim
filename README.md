# cophesim

Simulation tool to simulate phenotypes: continuous, dichotomous and survival for common variants from existing genotype data simulated with some other tool.

Version: 1.3.1

© Ilya Y. Zhbannikov, 2016-10-28

##Requirements

1. Python 2.7
2. plinkio https://pypi.python.org/pypi/plinkio
3. Plink v1.7 (to run examples)
4. R (to run examples)

## Installation

Download/clone the cophesim repository. The software is ready to use. You may add the home directory of cophesim to the PATH:

```
export PYTHONPATH=<cophesim home directory>:$PYTHONPATH
```

## Input

Files in one of the following data formats: Plink (.bed, .bim, .fam); ms, msms, msHot (plain text file), Genome (plain text file).

## Output:

* Phenotype file. This file is in text format and has the following suffices depending on the simulated phenotype trait: _pheno_bin.txt, _pheno_cont.txt, _pheno_surv.txt representing dichotomous (binary), quantitative (continuous) and survival phenotype.}
* Genotype file(s). Can be in the following formats: Plink (.bed, .bim, .fam). Other possible output format: blossoc (for BLOSSOC, suffices  .blossoc_pos, .blossoc_geno), qtdt (for QTDT, suffices .ped, .map, .dat), tassel (for Tassel, suffices .poly, .trait), emmax (for EMMAX, suffices  .emma_geno, .emma_pheno).}
* Summary statistics file. This is a plain text file which keeps the information about the run.}

## Usage
```
$python cophesim.py -i <path genotype files> -out <ouput prefix> [other parameters]
```

###Other parameters:

```
  -h, --help            show the help message and exit
  -i IDATA, --input IDATA
                        Path with prefix to your input file(s).
  -o OUTPUT_PREFIX, --output OUTPUT_PREFIX
                        Output prefix.
  -itype ITYPE          Input format: plink (for Plink, default), ms (for ms,
                        msms, msHot), genome (for Genome).
  -d, --dichotomous     A flag for dichotomous phenotype, True by default.
  -c, --continuous      A flag for continous phenotype, False by default.
  -s, --suvival         A flag to simulate survival phenotype, False by
                        default.
  -ce CEFF              A path to the file with effect of each causal SNP.
                        Must be in format: snp_index:effect. One snp per line
  -otype OTYPE          Indicates output format, default=plink. Other possible
                        output format: blossoc (for BLOSSOC), qtdt (for QTDT),
                        tassel (for Tassel), emmax (for EMMAX).
  -alpha ALPHA          An 'alpha' parameter for inverse probability equation
                        for the Gompertz hazard (see Bender at al., Generating
                        survival times to simulate Cox proportional hazards
                        models), 2005.
  -epi EPIFILE          File with interacting SNPs. One pair per line. Format:
                        snp1_index,snp2_index,effect
  -weib                 A flag to use Weibull distribution for survival
                        phenotype. True by default.
  -gomp                 A flag to use Gompertz distribution for survival
                        phenotype. False by default.
```

## Examples

```
plink --simulate-ncases 5000 --simulate-ncontrols 5000 --simulate wgas.sim --out sim.plink --make-bed

python cophesim.py -i sim.plink -o testout
```
The first command runs the data simulation. Here we simulate genetic dataset of 10k individuals, 5k cases and 5k controls. SNPs defined in \texttt{wgas.sim} (should be in the \texttt{cophesim} home directory). Then with next command we add a phenotype (dichotomous by default) to simulated genetic data.

To simulate continuous phenotypic trait, add the '-c' flag:

```
python cophesim.py -i sim.plink -o testout -c
```

This will simulate both continuous and dichotomous traits. To simulate survival trait, add '-s' flag:

```
python cophesim.py -i sim.plink -o testout -s
```

To test the output from other tools I prepared another repository: https://bitbucket.org/izhbannikov/cophesim_data . Clone or download this repository.

Example reading GENOME output:

```
python cophesim.py -i ../cophesim_data/tools/genome-0.2/genome.out.txt -o testout -s -c -itype genome -otype emmax
```

Example reading ms output:

```
python cophesim.py -i ../cophesim_data/tools/ms/ms.out.txt -o testout -s -c -itype ms -otype blossoc
```

