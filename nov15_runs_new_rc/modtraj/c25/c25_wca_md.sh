#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N cg_md
date

python md.py wca 1 1 25 1700 /home/cask0/home/tsanyal/c25ld/data/lammpstraj/c25/unconstrained_wca/c25_unbiased_mapped.lammpstrj.gz c25_wca_SP /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_c25/c25_wca_SP 7.8
python md.py wca 1 1 25 1700 /home/cask0/home/tsanyal/c25ld/data/lammpstraj/c25/unconstrained_wca/c25_unbiased_mapped.lammpstrj.gz c25_wca_SPLD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_c25/c25_wca_SPLD 7.8
python md.py wca 1 1 25 1700 /home/cask0/home/tsanyal/c25ld/data/lammpstraj/c25/unconstrained_wca/c25_unbiased_mapped.lammpstrj.gz c25_wca_LD /home/cask0/home/tsanyal/c25ld/nov15_runs_new_rc/cgff_c25/c25_wca_LD 7.8
