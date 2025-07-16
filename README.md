# gmxware

A collection of useful Gromacs mdp files that were most useful many years ago. Some modifications might be required to support newer versions of Gromacs.

The general protocol used for MD runs was:
1. Constrain the protein and execute
   - energy minimization
   - NVT simulation to gradually bring the solvent to room temperature
   - NPT simulation to let the solvent equilibrate

1. Release constraints and repeat the cycle without constraints.
   - energy minimization
   - NVT simulation to gradually bring the solvent to room temperature
   - NPT simulation to let the system equilibrate

1. Production run.

    Production run is just a very long 'NPT simulation to let the system equilibrate'.