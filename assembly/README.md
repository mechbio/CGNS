# assembly: assembly module

Date:               06-Jul-22
Version:            0.7.0

### Description of files in this directory
README.md           this file
param/*             calculations of model parameters
gen/*               generates initial configuration of nuclear components
sim/*               assembly of nuclear components
viz/*               visualization of the assembled structure
data/*              storage of numerical data
fig/*               save visualizations

### Usage
0. Please ignore the param/ directory for now.

1. The gen/README.md has instructions on generating the initial configuration
   (data/ne.data).

3. Next, the sim/README.md has instructions on simulating the nucleus.

4. Finally, the viz/README.md has instructions on visualization of the
   simulated nucleus.

5. The data/* directory (created at runtime) stores all numerical output. A
   description of files expected to be output in this directory is included
   below.

6. The visualizations can be saved to fig directory (created at runtime).

### Description of files in ./data/
ne.data             initially generated NE structure
init.out            LAMMPS log output during assembly simulation
init.lammpstrj      particle positions during assembly simulation
