#!/bin/bash
export OMP_NUM_THREADS=1
# for short simulations/debugging
time $PWD/../../Ndyn <$0 &>test.out &
# for long simulations
#nohup $PWD/gefs <$0 2>err >out &
exit
###################
START
#help		# will list directive options
###################
#Model 	ModelFCC #Model.xyz
BoxSize 15.0 15.0 15.0
NeiMax  70      # max number of neighbours
Types 1  #atom and particle types
 pt 1.d0
# Cu  63.546  #typeName and mass
Pairs 1
 pt pt sym cosine #test potential
OutFreq 40
ConfOutFreq 0 # configuration output frequency
Run 400 0.0001    # run for [i] steps each at [r] deltat
###################
END 
###################
