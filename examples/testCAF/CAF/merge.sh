#!/bin/bash

nn=`head -n 1 out.xyz_0001.dd`
m=`ls -1 out.xyz_0*.dd | wc -l | awk '{print $1}'`

n=`expr $nn \* $m`
echo $n
echo ""

for l in `ls out.xyz_0*.dd`; do
  tail -n $nn $l 
done
