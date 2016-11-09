#!/bin/bash
#-------------------------------------------Example begins----------------------------------------#
#Step 1: genetic data simulation:
plink --simulate-ncases 5000 --simulate-ncontrols 5000 --simulate wgas.sim --out sim.plink --make-bed
#Step 2: Convert .bed to .ped:
plink --bfile sim.plink --recode --out sim.plink
#Step3: phenotype simulation from previously made genetic data:
python cophesim.py -i sim.plink -o testout -itype plink -otype plink -c -ce effects.txt -s -gomp
