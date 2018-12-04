#!/bin/bash

#SBATCH -p shared
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH -J LRE_VELS
#SBATCH -C haswell

# Set directories
rundir=/global/homes/c/cschreck/production/radial_growth_damped_budding_basic
outdir=/global/cscratch1/sd/cschreck/radial_growth_damped_budding_basic

# geometric parameters
ar1=1.01
ar2=2.0

# rates
rate0=1d0
b=$1

# steps
steps=$2
layerskip=$3
dataskip=$4
prodskip=$5
restskip=$6
dt=$7

# growth layer widths
layerwidth=$8
layerdepth=$9
propdepth=${10}
bounddepth=${11}

# logical variables
movie=${12}
restart=${13}

# run parameters
desync=0.4
seed=-${14}

# output files
suffix=radial_layer${layerdepth}_desync${desync}_b${b}_seed${seed}_axial.dat
prodfile=prod_$suffix
restfile=restart_$suffix
radfile=radius_$suffix

cd $outdir

# run program
time $rundir/dimers_damped_radial_front_removecells_block_axial_spatialinfo_clean_v6.o <<EOF
  $ar1
  $ar2
  $rate0
  $b
  $steps
  $layerskip
  $dataskip
  $prodskip
  $restskip
  $dt
  $layerwidth
  $layerdepth
  $propdepth
  $bounddepth
  $desync
  $seed
  $prodfile 
  $restfile 
  $radfile 
  $movie
  $restart
EOF
