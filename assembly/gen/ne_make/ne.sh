#!/bin/bash
# -*- coding: utf-8 -*-

# This file is part of CGNS - Coarse-Grained Nucleus Simulator.
#
# Copyright 2022 Pranjal Singh
#
# When contributing, please append a new line (e.g. # Copyright [Year] [Name])
# to the above copyright notice.
#
# See the README file in the top-level CGNS directory.
# This software is distributed under the GNU General Public License.

# -----------------------------------------------------------------------------
# This file (ne.sh) is for generating init.data.
# -----------------------------------------------------------------------------

# date         :06-Jul-22
# version      :0.7.0
# usage        :./init.sh
# sh_version   :5.0.17(1)-release

# Run init.py
python3 ne.py 2

# copy to data
mkdir -p ../../data/
mv ne.data ../../data/

exit 0
