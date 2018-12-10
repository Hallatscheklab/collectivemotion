#!/bin/bash

# Set directories
outdir=~/collectivemotion_repo/cell_based_radial_budding/data/
rundir=~/collectivemotion_repo/cell_based_radial_budding/codes/

# geometric parameters
ar1=1.01
ar2=2.0

# rates
rate0=1d0
b=4e3

# steps
steps=1200000
layerskip=800
dataskip=5000
prodskip=5000
restskip=5000
dt=5e-6

# growth layer widths
layerwidth=2.1
layerdepth=9.0
propdepth=4.0
bounddepth=4.0

# logical variables
movie=.true.
restart=.true.

# run parameters
desync=0.4
seed=-1

# output files
suffix=radial_layer${layerdepth}_desync${desync}_b${b}_seed${seed}_axial.dat
prodfile=prod_$suffix
restfile=restart_$suffix
radfile=radius_$suffix

cd $outdir

# run program
time $rundir/dimers_damped_radial_removecells_axial.o <<EOF
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
