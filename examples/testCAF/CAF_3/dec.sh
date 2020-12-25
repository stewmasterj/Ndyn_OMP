#!/bin/bash

# decompose the model domain

f=model.xyz

p=/home/stewmasterj/progs/fortran/domainDecomposition

$p/domdec -N 3 3 3 -v -p -r 6.0 -b 65.322 65.322 65.322 ../$f

mv ../${f}_00*.dd .


# run this as a check 
echo "test 1"
cafrun -np 1 $p/testddCAF 6.0 T ../$f >test1.out
echo "test 27"
cafrun -np 27 $p/testddCAF 6.0 T $f >test27.out
