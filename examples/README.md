## Examples

### test
 randomly distributes 1000 particles into a 15x15x15 box. relaxes the configuration with a cos^2 repulsive potential.

### test\_larger
 larger case of *test*

### lammpsBench
 same as *test\_larger* but using LAMMPS as a benchmark

### SineTest
 uses a soft potential for relaxation reading a LAMMPS dump file as input.

### SineTestLMP
 same as *SineTest* but using LAMMPS as a benchmark

### MorseTest
 reads a copper FCC crystal from LAMMPS dump file
 
### MorseTestLmp
 same as *MorseTest* but using LAMMPS as a benchmark

### testPDspringTension
 uses springs based on interparticle strain for force calculation, like Peridynamics.
 and has a critical strain criterion.

### testFollow
 This test should be the same as *testPDspringTension* except that it has an
 extra particle that simply follows the deformation field.

### testT
 uses the Copper FCC model to test various ensembles and tension loads.



