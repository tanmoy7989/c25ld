#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N methane_clust
date

#wca with 7.8 cutoff for AA trajectories
python methane_cluster_analyse.py ~/c25ld/data/lammpstraj/methane_wca_transferability_runs/methane10_wca_mapped.lammpstrj.gz wca 10 methane10_wca_highrc_AA geom_prop/methane_LDKnot/methane10
python methane_cluster_analyse.py ~/c25ld/data/lammpstraj/methane_wca_transferability_runs/methane20_wca_mapped.lammpstrj.gz wca 20 methane20_wca_highrc_AA geom_prop/methane_LDKnot/methane20
python methane_cluster_analyse.py ~/c25ld/data/lammpstraj/methane_waterbox/methane30_wca_mapped.lammpstrj.gz wca 30 methane30_wca_highrc_AA geom_prop/methane_LDKnot/methane30


# wca with 7.8 cutoff for CG trajectories
python methane_cluster_analyse.py modtraj/methane_LDKnot/methane25/methane25_wca_highrc_SP.lammpstrj.gz wca 25 methane25_wca_highrc_SP geom_prop/methane_LDKnot/methane25
python methane_cluster_analyse.py modtraj/methane_LDKnot/methane25/methane25_wca_highrc_SPLD.lammpstrj.gz wca 25 methane25_wca_highrc_SPLD geom_prop/methane_LDKnot/methane25
python methane_cluster_analyse.py modtraj/methane_LDKnot/methane25/methane25_wca_highrc_LD.lammpstrj.gz wca 25 methane25_wca_highrc_LD geom_prop/methane_LDKnot/methane25

# wca transferability with 7.8 cutoff for CG trajectories
python methane_cluster_analyse.py modtraj/methane_LDKnot/methane10/methane10_wca_highrc_SP.lammpstrj.gz wca 10 methane10_wca_highrc_SP geom_prop/methane_LDKnot/methane10
python methane_cluster_analyse.py modtraj/methane_LDKnot/methane10/methane10_wca_highrc_SPLD.lammpstrj.gz wca 10 methane10_wca_highrc_SPLD geom_prop/methane_LDKnot/methane10
python methane_cluster_analyse.py modtraj/methane_LDKnot/methane10/methane10_wca_highrc_LD.lammpstrj.gz wca 10 methane10_wca_highrc_LD geom_prop/methane_LDKnot/methane10

python methane_cluster_analyse.py modtraj/methane_LDKnot/methane20/methane20_wca_highrc_SP.lammpstrj.gz wca 20 methane20_wca_highrc_SP geom_prop/methane_LDKnot/methane20
python methane_cluster_analyse.py modtraj/methane_LDKnot/methane20/methane20_wca_highrc_SPLD.lammpstrj.gz wca 20 methane20_wca_highrc_SPLD geom_prop/methane_LDKnot/methane20
python methane_cluster_analyse.py modtraj/methane_LDKnot/methane20/methane20_wca_highrc_LD.lammpstrj.gz wca 20 methane20_wca_highrc_LD geom_prop/methane_LDKnot/methane20

python methane_cluster_analyse.py modtraj/methane_LDKnot/methane30/methane30_wca_highrc_SP.lammpstrj.gz wca 30 methane30_wca_highrc_SP geom_prop/methane_LDKnot/methane30
python methane_cluster_analyse.py modtraj/methane_LDKnot/methane30/methane30_wca_highrc_SPLD.lammpstrj.gz wca 30 methane30_wca_highrc_SPLD geom_prop/methane_LDKnot/methane30
python methane_cluster_analyse.py modtraj/methane_LDKnot/methane30/methane30_wca_highrc_LD.lammpstrj.gz wca 30 methane30_wca_highrc_LD geom_prop/methane_LDKnot/methane30


# lj with 6.5 cutoff
python methane_cluster_analyse.py ~/c25ld/data/lammpstraj/methane25/nvt/methane_lj_unbiased_mapped.lammpstrj.gz lj 25 methane25_lj_AA geom_prop/methane25
python methane_cluster_analyse.py ../data/modtraj/old_runs/methane/unconstrained_lj/SP/methane_lj_SP.lammpstrj.gz lj 25 methane25_lj_SP geom_prop/methane25
python methane_cluster_analyse.py ../data/modtraj/old_runs/methane/unconstrained_lj/SPLD/methane_lj_SPLD.lammpstrj.gz lj 25 methane25_lj_SPLD geom_prop/methane25
python methane_cluster_analyse.py ../data/modtraj/old_runs/methane/unconstrained_lj/LD/methane_lj_LD.lammpstrj.gz lj 25 methane25_lj_LD geom_prop/methane25
