#!/bin/bash
# -*- coding: utf-8 -*-

# CGNS - Coarse-Grained Nucleus Simulator
#
# See the README file in the top-level CGNS directory.
# This software is released under the GNU General Public License.

# -----------------------------------------------------------------------------
# This file (ne.sh) is for generating init.data.
# -----------------------------------------------------------------------------

# date         :06-Jul-22
# version      :0.7.0
# usage        :./init.sh
# sh_version   :5.0.17(1)-release

# Run init.py
python3 ne.py 1.5

# copy to data
mkdir -p ../../data/
mv ne.data ../../data/

exit 0
