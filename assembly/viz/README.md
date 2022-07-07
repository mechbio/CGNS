# viz: visualization sub-module

Date:               08-Jul-22
Version:            0.7.0

### Description of files in this directory
README.md             this file
viz.ovito             visualization of ne.data and init.lammpstrj in ../data/
viz.sh                processing structure file for visualization (viz)
composition.ovito     viz of CGNP (Fig. 1a) and its composition (Fig. 1b)
slice-diagonal.ovito  viz of a slice of CGNP along the xy-diagonal (Fig. 1c)
surface.ovito         viz of nuclear surface deformation (Fig. 2a)

### Usage
0. Only the viz.ovito file is operational currently. Please use only this file.
1. On first time run, the OVITO file may require manual reloading of
   ne.data and init.lammpstrj from ../data/ into the 'Data source' and the
   'Load trajectory' in the 'Modifications' pipeline respectively.
