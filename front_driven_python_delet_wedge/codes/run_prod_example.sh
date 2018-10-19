#!/bin/bash
 
# input parameters
T=5.0
L=600.0
w0=240.0
diam=1.0
ndis=4
diff=0.0
steps=600000
skipb=500
skiph=5000
dh0=0.001
s=-0.06
seed=1  
trials=1

# output file names
suff=L${L}_T${T}_s${s}_diff${diff}_seed${seed}.dat
fileb=bound_$suff
fileh=front_$suff

python front_driven_prod.py <<EOF
  $T
  $L
  $w0
  $diam
  $ndis
  $diff
  $steps
  $skipb
  $skiph
  $dh0
  $s
  $seed
  $trials
  '$fileb'
  '$fileh'
EOF
