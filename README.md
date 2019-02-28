#  N D Y N
FORTRAN90/OpenMP version 0.10.18

*A newtonian dynamics program.*

written by Ross J. Stewart


This is a molecular dynamic program that is multithreaded for use on shared memory
 machines. This does not use domain decomposition like most MPI MD codes do.
The various threads simply share the primary particle loops.
This code only utilizes particle pair potentials, no multibody potentials have
 been enabled.

# Preferred units 
length:         1 Angstrom\
energy:         1 electron volt [eV]\
mass:           1 atomic mass unit [amu]\
time:           10.1805057e-15 (sec) 10.1805057 [fs]\
temperature:    1 Kelvin [K]\
pressure:       1 eV.Ang^-3 = 160.217662080e9 Pa\
Boltzmann Constant: k=8.6173303e-5 eV/K\
charge:         1 e.3.794685 = e.sqrt(14.399637 Ang.eV/e^2) : [e.sqrt(Ang.eV)]\

To convert total KE to T\
T = 2KE/(3Nk)\
KE/N = 3Tk/2
