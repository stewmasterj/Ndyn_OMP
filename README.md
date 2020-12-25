#  N D Y N
FORTRAN90/OpenMP         version 0.10.18  
FORTRAN90/OpenMP-Coarray version 0.12.20  

*A newtonian dynamics program.*

written by Ross J. Stewart


This is a molecular dynamic program that is multithreaded for use on shared memory
 machines. 
There are two versions of this program, the OpenMP only version in ./src/ and
 a hybrid parallel version using both OpenMP for shared memory and Open-Coarray-Fortran
 for distributed memory on massively parallel archetectures in ./src\_CAF.
Open-Coarray-Fortran (CAF) is in some ways a compiler-level wrapper around the MPI library.
The CAF version requires the model domain to be decomposed using the code in ../DomDec.
Many of the included examples makes use of the crystal tesselation coce in ../tessel.
This code only utilizes particle pair potentials, no multibody potentials like Tersoff 
 or EAM have been enabled.

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
