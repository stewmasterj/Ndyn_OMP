#!/bin/bash

# decompose the model domain

f=model.xyz

p=/home/stewmasterj/progs/fortran/domainDecomposition

$p/domdec -N 2 1 1 -v -p -r 6.0 -b 65.322 65.322 65.322 ../$f

mv ../${f}_00*.dd .


# run this as a check 
echo "test 1"
cafrun -np 1 $p/testddCAF 6.0 T ../$f >test1.out
echo "test 2"
cafrun -np 2 $p/testddCAF 6.0 T $f >test2.out
