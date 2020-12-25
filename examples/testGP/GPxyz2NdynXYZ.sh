#!/bin/bash

# convert the XYZ files that mkXtal generates and place the scale number as
#   part of the atom type.

if [ -z $1 ]; then
   echo "please specify a file name to convert"
   exit
fi

N=`head  -n 1 $1 | awk '{print $1}'`

#writting new output to stdout
echo $N
echo ""
tail -n+3 $1 | awk '{print $1$7"\t"$2"\t"$3"\t"$4}'


