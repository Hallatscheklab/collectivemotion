#!/bin/bash

dir=/global/homes/c/cschreck/production/radial_growth_damped_budding_basic

#att=0.0

b=4e3
#steps=1200000
steps=1270000
layerskip=800
dataskip=5000
prodskip=5000
restskip=5000
#treeskip=5000
dt=5e-6

layerwidth=2.1
layerdepth=9.0
propdepth=4.0
bounddepth=4.0

movie=.true.
restart=.true.

seed=1
#$dir/run_radial_v5.sh $att $b $steps $layerskip $dataskip $prodskip $restskip $treeskip $dt $layerwidth $layerdepth $propdepth $bounddepth $movie $restart $seed
$dir/run_radial_v6.sh $b $steps $layerskip $dataskip $prodskip $restskip $dt $layerwidth $layerdepth $propdepth $bounddepth $movie $restart $seed
