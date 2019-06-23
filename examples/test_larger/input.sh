#!/bin/bash
export OMP_NUM_THREADS=1
# for short simulations/debugging
time $PWD/../../Ndyn <$0 &>out &
exit
###################
START
#help		# will list directive options
###################
BoxSize 40.0 40.0 40.0 #roughly 20 particles per length
NTotal 64000  #number of particles in random box
#Model 	ModelFCC #Model.xyz
NeiMax  210      # max number of neighbours
Types 1  #atom and particle types
 pt 1.d0
Pairs 1
 pt pt sym cosine #test potential
OutFreq 40
ConfOutFreq 0 # configuration output frequency
Run 4000 0.0001    # run for [i] steps each at [r] deltat
###################
END 
###################
