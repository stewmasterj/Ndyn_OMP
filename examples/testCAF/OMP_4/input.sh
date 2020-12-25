#!/bin/bash
export OMP_NUM_THREADS=4
# for short simulations/debugging
time $PWD/../../../Ndyn <$0 &>out &
exit
###################
START
#help		# will list directive options
###################
BoxSize 65.322 65.322 65.322 
NeiMax  200      # max number of neighbours
Types 1  #atom and particle types
 Cu 63.546 
Skin 0.0 
Pairs 1
 Cu Cu sym morse 6.0 0.3429 1.359 2.866
Model 	../model.xyz xyz
#WritePotentialTable
OutFreq 10 #40
ConfOutFreq 100 xyz # configuration output frequency
Run 100 0.0001    # run for [i] steps each at [r] deltat
###################
END 
###################
