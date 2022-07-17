# viz: visualization sub-module

Date:               17-Jul-22
Version:            0.9.0

### Description of files in this directory
README.md             this file
viz.sh                processing structure file for visualization
viz3.ovito            visualize NE vesicle made via ../sim/init3.lmp
viz2.ovito            visualize NM vesicle made via ../sim/init2.lmp
composition.ovito     from CGNPS; kept for later use
slice-diagonal.ovito  from CGNPS; kept for later use
surface.ovito         from CGNPS; kept for later use

### Usage
0. Either of viz2 and viz3 OVITO files can be used.
1. On first time run, the OVITO file may require manual reloading of
   ne.data and init.lammpstrj from ../data/ into the 'Data source' and the
   'Load trajectory' in the 'Modifications' pipeline respectively.
