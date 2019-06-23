#!/bin/bash
export OMP_NUM_THREADS=2
# for short simulations/debugging
$PWD/../../Ndyn <$0 &>out &
# for long simulations
#nohup $PWD/gefs <$0 2>err >out &
exit
###################
START
#help		# will list directive options
###################
#BoxSize 40.0 40.0 40.0 #roughly 20 particles per length
#NTotal 64000  #number of particles in random box
Model 	model.dump dump
#Cutoff	6.0  #cutoff radius
Seed	1234567890
Skin 	0.0
#NeiMax  35      # max number of neighbours
Types 2  #atom and particle types
 Cu   63.546 0.0 # name mass charge
 Cu-1 63.546 0.0 # name mass charge
#WritePotentialTable
Pairs 2
 Cu Cu sym spring 6.0  2.152295e-3 0.1 12.459008 #rc c s0 Vi
 Cu Cu-1 asym spring 6.0  2.152295e-3 0.1 12.459008 #rc c s0 Vi
# Cu Cu sym morse 6.0  0.3429 1.359 2.866
# stiffness is c= 12E/(pi*rc^4)
# strain    is s0= 1/6*sqrt(30*G/(E*rc))
# then  dt  is dt= 0.8*sqrt(2*dens/(pi*rc^2*bc))

# run 1 NP to minimize configurational energy
OutFreq 	100
ReNeighbor 	F
Exec 		rm out.dump
ConfOutFreq 	100 dump # configuration output frequency
BoxVel		1.e-2 0.0 0.0 #load the box at this rate Ang/(10 fs)
Run		40000 0.1
###################
END 
###################
