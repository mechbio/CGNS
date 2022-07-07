#!/bin/bash
# -*- coding: utf-8 -*-

# CGNS - Coarse-Grained Nucleus Simulator
#
# See the README file in the top-level CGNS directory.
# This software is released under the GNU General Public License.

# -----------------------------------------------------------------------------
# This file (run.sh) is used to run and monitor simulations.
# -----------------------------------------------------------------------------

# date         :06-Jul-22
# version      :0.7.0
# usage        :./run.sh
# sh_version   :5.0.17(1)-release

# Ensure output file exists
touch ../data/init.out
# Parallel run
mpirun -np 4 lmp_mylammps < init2.lmp > ../data/init.out
# Status
# watch -n 1 "tail -10 ../data/init.out"

exit 0
