#!/bin/bash
export OMP_NUM_THREADS=1
# for short simulations/debugging
$PWD/../../../Ndyn <$0 &>out &
# for long simulations
#nohup $PWD/gefs <$0 2>err >out &
exit
###################
START
#help		# will list directive options
###################
BoxSize 87.096 43.548 43.548 # 12*2*a0, 4*2*a0, 4*2*a0
#Model 	../MorseTest/model.dump dump
Seed	1234567890
Skin 	0.2
#NeiMax  300      # max number of neighbours
Types 1 #4  #atom and particle types
 Cu2   1016.736 0.0 # name mass charge
# Cu1   63.546 0.0 # name mass charge
# Cu-1  63.546 0.0 # name mass charge
# Cu-2 1016.736 0.0 # name mass charge
Model	model.xyz  #named atom types must read model after type definitions
#WritePotentialTable
Pairs 1 #6
 Cu2  Cu2   sym Morse  12.0  5.4864 0.6795 5.732  # 16xD0 for proper stiffness
 #Cu2  Cu2   sym Morse  12.0  2.7432 0.96096 5.732  # 8xD0 alpha/sqrt(2)
 #Cu2  Cu2   sym Morse  12.0  2.7432 0.6795 5.732  # 8xD0
# Cu1  Cu1   sym Morse  6.0  0.3429 1.359 2.866
# Cu-1 Cu1  asym Morse  6.0  0.3429 1.359 2.866
# Cu1  Cu-2 asym spring 6.0  2.152295e-3 0.1 12.459008 #rc c s0 Vi
# Cu-2 Cu2  asym Morse  12.0  2.7432 0.6795 5.732
# Cu2  Cu-1 asym spring 12.0  2.152295e-3 0.1 12.459008 #rc c s0 Vi
# Cu Cu sym spring 6.0  2.152295e-3 0.1 12.459008 #rc c s0 Vi
# Cu Cu sym morse 6.0  0.3429 1.359 2.866
# stiffness is c= 12E/(pi*rc^4)
# strain    is s0= 1/6*sqrt(30*G/(E*rc))
# then  dt  is dt= 0.8*sqrt(2*dens/(pi*rc^2*bc))

# output options
OutFreq 	10
ReNeighbor 	F # shouldn't need to reneighbour, saves a little time checking
Exec 		rm out.dump
ConfOutFreq 	10 dump # configuration output frequency

### Run parameters ###
# run 1 relax stresses with a barostat at low T
Isobaric        0.0 0.0 0.0  0.1 10
Run             1000 0.1
Anisobaric      #unset the barostat (constant volume)
# run 2 NVT with initial velocity distribution
Velocity        0.3102238908  #3/2*k*300K, k=8.6173303e-5 eV/K
Isokinetic      0.3102238908 0.1 10
ReNeighbor      T #might need this if skin is too small
Run             1000 0.1    # run for [i] steps each at [r] deltat
# run 3 NPT
Isobaric        0.0 0.0 0.0  1.0 10
Run             2000 0.1
Anisokinetic    #unset thermostat (constant E)
Anisobaric      #unset the barostat (constant volume)
# stretch the box
ReNeighbor 	T #F #stretching will need reneighbouring if any reorientations
BoxVel		1.e-2 0.0 0.0 #load the box at this rate Ang/(10 fs)
Run		5000 0.1
###################
END 
###################
