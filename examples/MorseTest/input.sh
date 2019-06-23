#!/bin/bash
export OMP_NUM_THREADS=3
# for short simulations/debugging
time $PWD/../../Ndyn <$0 &>out &
exit
###################
START
#help		# will list directive options
###################
Model 	model.dump dump
NeiMax  200      # max number of neighbours
Types 1  #atom and particle types
 Cu 63.546 
Pairs 1
 Cu Cu sym morse 6.0 0.3429 1.359 2.866
#WritePotentialTable
OutFreq 40
ConfOutFreq 0 dump # configuration output frequency
Run 4000 0.0001    # run for [i] steps each at [r] deltat
###################
END 
###################
