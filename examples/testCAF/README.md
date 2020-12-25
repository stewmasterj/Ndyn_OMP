first build the crystal we will perform tests on.  

We want to compare the Ndyn\_OMP with Ndyn\_CAF to make sure they're the same.  
So we need to compare using a large crystal model, able to decompose.  

1) run the mkXtal script named 'in.mX'  
   Be sure that the unitcell file is correctly located in the in.mX file.  
      mkXtal < in.mX | tee mX.out  

2) run the Ndyn\_OMP code on the model in DIR=OMP\_4  

3) decompose the model domain using 'dec.sh' in each subdirectory  
  3a) 2x1x1  OMP=2 DIR=CAF  
  3b) 2x2x2  OMP=1 DIR=CAF\_2  
  3c) 3x3x3  OMP=1 DIR=CAF\_3  

4) run the cases with Ndyn\_CAF using 'input.sh'  
  the results in 'out' should basically be the same.  


NOTES:  
cannot have SYNC ALL statements inside OMP parallel constructs
