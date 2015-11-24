#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N geom_analyse

date
################ c25 #########################
python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/c25/unconstrained_lj/c25_unbiased_mapped.lammpstrj.gz "(lj,AA)" 25 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c25

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/c25/unconstrained_wca/c25_unbiased_mapped.lammpstrj.gz "(wca,AA)" 25 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c25

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/c25/unconstrained_lj/SP/c25_lj_SP.lammpstrj.gz "(lj,SP)" 25 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c25

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/c25/unconstrained_lj/SPLD/c25_lj_SPLD.lammpstrj.gz "(lj,SPLD)" 25 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c25

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/c25/unconstrained_lj/LD/c25_lj_LD.lammpstrj.gz "(lj,LD)" 25 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c25

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/c25/unconstrained_wca/SP/c25_wca_SP.lammpstrj.gz "(wca,SP)" 25 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c25

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/c25/unconstrained_wca/SPLD/c25_wca_SPLD.lammpstrj.gz "(wca,SPLD)" 25 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c25

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/c25/unconstrained_wca/LD/c25_wca_LD.lammpstrj.gz "(wca,LD)" 25 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c25




date
################ c10 #########################
python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/polymer_transferability_runs/c10/unconstrained_lj/c10_unbiased_mapped.lammpstrj.gz "(lj,AA)" 10 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c10

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/polymer_transferability_runs/c10/unconstrained_wca/c10_unbiased_mapped.lammpstrj.gz "(wca,AA)" 10 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c10

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c10/unconstrained_lj/SP/c10_lj_SP.lammpstrj.gz "(lj,SP)" 10 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c10

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/aug15_runs_c10/unconstrained_lj/SPLD/c10_lj_SPLD.lammpstrj.gz "(lj,SPLD)" 10 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c10

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c10/unconstrained_lj/LD/c10_lj_LD.lammpstrj.gz "(lj,LD)" 10 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c10

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c10/unconstrained_wca/SP/c10_wca_SP.lammpstrj.gz "(wca,SP)" 10 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c10

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c10/unconstrained_wca/SPLD/c10_wca_SPLD.lammpstrj.gz "(wca,SPLD)" 10 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c10

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c10/unconstrained_wca/LD/c10_wca_LD.lammpstrj.gz "(wca,LD)" 10 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c10




date
################ c20 #########################
python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/polymer_transferability_runs/c20/unconstrained_lj/c20_unbiased_mapped.lammpstrj.gz "(lj,AA)" 20 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c20

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/polymer_transferability_runs/c20/unconstrained_wca/c20_unbiased_mapped.lammpstrj.gz "(wca,AA)" 20 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c20

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c20/unconstrained_lj/SP/c20_lj_SP.lammpstrj.gz "(lj,SP)" 20 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c20

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c20/unconstrained_lj/SPLD/c20_lj_SPLD.lammpstrj.gz "(lj,SPLD)" 20 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c20

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c20/unconstrained_lj/LD/c20_lj_LD.lammpstrj.gz "(lj,LD)" 20 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c20

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c20/unconstrained_wca/SP/c20_wca_SP.lammpstrj.gz "(wca,SP)" 20 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c20

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c20/unconstrained_wca/SPLD/c20_wca_SPLD.lammpstrj.gz "(wca,SPLD)" 20 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c20

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c20/unconstrained_wca/LD/c20_wca_LD.lammpstrj.gz "(wca,LD)" 20 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c20




date
################ c30 #########################
python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/polymer_transferability_runs/c30/unconstrained_lj/c30_unbiased_mapped.lammpstrj.gz "(lj,AA)" 30 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c30

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/polymer_transferability_runs/c30/unconstrained_wca/c30_unbiased_mapped.lammpstrj.gz "(wca,AA)" 30 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c30

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c30/unconstrained_lj/SP/c30_lj_SP.lammpstrj.gz "(lj,SP)" 30 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c30

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c30/unconstrained_lj/SPLD/c30_lj_SPLD.lammpstrj.gz "(lj,SPLD)" 30 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c30

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c30/unconstrained_lj/LD/c30_lj_LD.lammpstrj.gz "(lj,LD)" 30 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c30

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c30/unconstrained_wca/SP/c30_wca_SP.lammpstrj.gz "(wca,SP)" 30 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c30

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c30/unconstrained_wca/SPLD/c30_wca_SPLD.lammpstrj.gz "(wca,SPLD)" 30 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c30

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c30/unconstrained_wca/LD/c30_wca_LD.lammpstrj.gz "(wca,LD)" 30 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c30




date
################ c40 #########################
python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/polymer_transferability_runs/c40/unconstrained_lj/c40_unbiased_mapped.lammpstrj.gz "(lj,AA)" 40 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c40

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/polymer_transferability_runs/c40/unconstrained_wca/c40_unbiased_mapped.lammpstrj.gz "(wca,AA)" 40 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c40

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c40/unconstrained_lj/SP/c40_lj_SP.lammpstrj.gz "(lj,SP)" 40 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c40

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c40/unconstrained_lj/SPLD/c40_lj_SPLD.lammpstrj.gz "(lj,SPLD)" 40 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c40

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c40/unconstrained_lj/LD/c40_lj_LD.lammpstrj.gz "(lj,LD)" 40 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c40

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c40/unconstrained_wca/SP/c40_wca_SP.lammpstrj.gz "(wca,SP)" 40 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c40

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c40/unconstrained_wca/SPLD/c40_wca_SPLD.lammpstrj.gz "(wca,SPLD)" 40 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c40

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c40/unconstrained_wca/LD/c40_wca_LD.lammpstrj.gz "(wca,LD)" 40 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c40




date
################ c50 #########################
python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/polymer_transferability_runs/c50/unconstrained_lj/c50_unbiased_mapped.lammpstrj.gz "(lj,AA)" 50 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c50

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/lammpstraj/polymer_transferability_runs/c50/unconstrained_wca/c50_unbiased_mapped.lammpstrj.gz "(wca,AA)" 50 1 1700 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c50

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c50/unconstrained_lj/SP/c50_lj_SP.lammpstrj.gz "(lj,SP)" 50 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c50

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c50/unconstrained_lj/SPLD/c50_lj_SPLD.lammpstrj.gz "(lj,SPLD)" 50 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c50

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c50/unconstrained_lj/LD/c50_lj_LD.lammpstrj.gz "(lj,LD)" 50 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c50

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c50/unconstrained_wca/SP/c50_wca_SP.lammpstrj.gz "(wca,SP)" 50 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c50

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c50/unconstrained_wca/SPLD/c50_wca_SPLD.lammpstrj.gz "(wca,SPLD)" 50 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c50

python geom_analyse.py /home/cask0/home/tsanyal/c25ld/data/modtraj/polymer_transferability_runs/c50/unconstrained_wca/LD/c50_wca_LD.lammpstrj.gz "(wca,LD)" 50 1 0000 /home/cask0/home/tsanyal/c25ld/data/analysis/geom_prop/c50
