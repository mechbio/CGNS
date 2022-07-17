#!/bin/bash
# -*- coding: utf-8 -*-

# This file is part of CGNES - Coarse-Grained Nuclear Envelope Simulator.
#
# Copyright 2021, 2022 Pranjal Singh
#
# When contributing, please append a new line (e.g. # Copyright [Year] [Name])
# to the above copyright notice.
#
# See the README file in the top-level CGNES directory.
# This software is distributed under the GNU General Public License.

# -----------------------------------------------------------------------------
# This file (viz.sh) is used to process data for visualization with OVITO.
# -----------------------------------------------------------------------------

# date         :09-Jul-22
# version      :0.9.0
# usage        :./viz.sh
# sh_version   :5.0.17(1)-release

# Make fig directory if it doesn't already exist
mkdir -p ../fig/

# Go to data
cd ../data

# Remove portion that is not handled in ovito-basic-3.3.5
sed -n '/Ellipsoids/q;p' < ne.data > ne2.data
# sed -n '/Ellipsoids/q;p' < init.data > init2.data

# Return to pwd
cd ../viz

exit 0
