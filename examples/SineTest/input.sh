#!/bin/bash
export OMP_NUM_THREADS=1
# for short simulations/debugging
time $PWD/../../Ndyn <$0 &>out &
exit
###################
START
#help		# will list directive options
###################
Model 	model.dump dump
#NeiMax  35      # max number of neighbours
Skin	1.0
Types 1  #atom and particle types
 Cu 63.546 
Pairs 1
 Cu Cu sym soft 6.0 0.5 6.0 
# Cu Cu sym mors 0.3429 1.359 2.866
#WritePotentialTable
OutFreq 40
ConfOutFreq 0 dump # configuration output frequency
Run 4000 0.00001    # run for [i] steps each at [r] deltat
###################
END 
###################
