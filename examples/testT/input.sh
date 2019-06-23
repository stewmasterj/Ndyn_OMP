#!/bin/bash
export OMP_NUM_THREADS=2
# for short simulations/debugging
$PWD/../../Ndyn <$0 &>out &
exit
###################
START
###################
Model 	../MorseTest/model.dump dump
Seed	1234567890
Skin 	0.2
#NeiMax  35      # max number of neighbours
Types 1  #atom and particle types
 Cu 63.546 0.0 # name mass charge
WritePotentialTable
Pairs 1
 Cu Cu sym morse 6.0  0.3429 1.359 2.866

# run 1 NP to minimize configurational energy
OutFreq 	100
ReNeighbor 	F
Exec 		rm out.dump
ConfOutFreq 	100 dump # configuration output frequency
Isobaric 	0.0 0.0 0.0  0.1 10
Run 		2000 0.1
Anisobaric 	#unset the barostat (constant volume)

# run 2 NVT with initial velocity distribution
Velocity 	0.03877798635  #3/2*k*300K, k=8.6173303e-5 eV/K
Isokinetic 	0.03877798635 0.1 10
ReNeighbor 	T
Run 		4000 0.1    # run for [i] steps each at [r] deltat

# run 3 NPT
Isobaric        0.0 0.0 0.0  0.1 10
Run		4000 0.1
Anisokinetic 	#unset thermostat (constant E)
Anisobaric 	#unset the barostat (constant volume)

# run 4 NVE uncontrolled
Run 		10000 0.1    # run for [i] steps each at [r] deltat

# run 5 NVE uncontrolled
BoxVel		0.0 0.004 0.0 #load the box in Y at 0.1 Ang/(10 fs)
Run		40000 0.2
###################
END 
###################
