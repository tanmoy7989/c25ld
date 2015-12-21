#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N methane_clust
date

# first for wca with 6.5 cutoff
#python methane_cluster_analyse.py ~/c25ld/data/lammpstraj/methane25/nvt/methane_wca_unbiased_mapped.lammpstrj.gz wca 25 methane25_wca_AA geom_prop/methane25
#python methane_cluster_analyse.py modtraj/methane25/methane25_wca_SP.lammpstrj.gz wca 25 methane25_wca_SP geom_prop/methane25
#python methane_cluster_analyse.py modtraj/methane25/methane25_wca_SPLD.lammpstrj.gz wca 25 methane25_wca_SPLD geom_prop/methane25
#python methane_cluster_analyse.py modtraj/methane25/methane25_wca_LD.lammpstrj.gz wca 25 methane25_wca_LD geom_prop/methane25

# then lj with 6.5 cutoff
python methane_cluster_analyse.py ~/c25ld/data/lammpstraj/methane25/nvt/methane_lj_unbiased_mapped.lammpstrj.gz lj 25 methane25_wca_lj geom_prop/methane25
python methane_cluster_analyse.py modtraj/methane25/methane25_lj_SP.lammpstrj.gz lj 25 methane25_lj_SP geom_prop/methane25
python methane_cluster_analyse.py modtraj/methane25/methane25_lj_SPLD.lammpstrj.gz lj 25 methane25_lj_SPLD geom_prop/methane25
python methane_cluster_analyse.py modtraj/methane25/methane25_lj_LD.lammpstrj.gz lj 25 methane25_lj_LD geom_prop/methane25
