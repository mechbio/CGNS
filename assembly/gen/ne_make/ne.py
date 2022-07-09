#!/usr/bin/env python3
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
# This file (ne.py) generates configuration data file (ne.data)
# which will be read by init.lmp.
# -----------------------------------------------------------------------------

# date         :06-Jul-22
# version      :0.7.0
# usage        :python3 ne.py 2.5
# py_version   :3.8.10

# Initialize ==================================================================

import sys
import numpy as np
import random

scl = float(sys.argv[1]) # scale of simulation box
L = 100*scl # Length of sides of simulation box in x and y (Lx = Ly = L)
sd = 221 # seed for random generator

m = 1.237 # mass of vesicle particles in reduced units

sclsq = scl*scl

# Simulation box
xlo = -L/2.0
xhi = L/2.0
ylo = -L/2.0
yhi = L/2.0
zlo = -L/2.0
zhi = L/2.0

# Initialize counts
natoms = 0
nbonds = 0
nangles = 0
nellipsoids = 0
natomTs = 0
nbondTs = 0
nangleTs = 0
nmols = 0

# Functions ===================================================================

def spherical_monolayer(D,d):
    """A spherical surface of diameter D with points spaced by d.
    Based on the points_on_sphere routine in fluidmembrane package.
    """

    R1 = D/2
    xyz = np.asarray([[0.0, 0.0, R1],
                        [0.0, 0.0, -R1]])

    for i in range(1,int(np.pi/2.0/np.arcsin(d0/2.0/R1))):
        ti = (i-1)*np.pi/int(np.pi/2.0/np.arcsin(d0/2.0/R1))
        for j in range(0,int(np.pi*np.sin(ti)/np.arcsin(d0/2.0/R1))):
            fj = (j-1)*2.0*np.pi/int(np.pi*np.sin(ti)/np.arcsin(d0/2.0/R1))
            x0 = R1*np.sin(ti)*np.cos(fj)
            y0 = R1*np.sin(ti)*np.sin(fj)
            z0 = R1*np.cos(ti)
            xyz = np.vstack([xyz, [x0, y0, z0]])

    quat = np.asarray([[0,0,0,0]])
    natoms = np.shape(xyz)[0]
    nxx = [1, 0, 0]
    for j in range(natoms):
        n = xyz[j,:]/np.linalg.norm(xyz[j,:]);
        nv = np.cross(nxx,n);
        nquat = nv/np.linalg.norm(nv) ;
        thetaquat=np.arccos(np.dot(nxx,n));
        quat = np.vstack([quat, [np.cos(thetaquat/2.0),nquat[0]*np.sin(thetaquat/2.0),nquat[1]*np.sin(thetaquat/2.0),nquat[2]*np.sin(thetaquat/2.0)]])

    return natoms, np.delete(quat, 0, axis=0), xyz

# Molecules ===================================================================

# spherical_monolayer ---------------------------------------------------------
d0 = 2.0
d1 = d0*1.08

natoms0, quat0, atoms0 = spherical_monolayer(0.8*L,d1)
# print(np.shape(atoms0))
# print(np.shape(quat0))
# data = atoms0
# x = data[:, 0]
# y = data[:, 1]
# z = data[:, 2]
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import proj3d
# fig = plt.figure(figsize=(8, 8))
# ax = fig.add_subplot(111, projection='3d')
# ax.scatter(x, y, z)
# plt.show()

types0 = 1
mass = m
shape0 = [d0, d0, d0]
volm = (4.0/3.0*(np.pi*shape0[0]*shape0[1]*shape0[2])/8.0)
density0 = mass/volm
ellipsoidflag0 = 1
molecule_ID0 = 0
charge0 = 0

natoms += natoms0
nellipsoids += natoms0
natomTs += 1

print(f'{natoms} spherical_monolayer particles generated.')

# Write =======================================================================

with open('ne.data','w') as fout:
    # Title comment -----------------------------------------------------------
    fout.write('Initial configuration of spherical_monolayer\n\n')

    # Header ------------------------------------------------------------------
    fout.write(f'{natoms} atoms\n')
    fout.write(f'{natomTs} atom types\n')
    fout.write(f'{nbonds} bonds\n')
    fout.write(f'{nbondTs} bond types\n')
    fout.write(f'{nangles} angles\n')
    fout.write(f'{nangleTs} angle types\n')
    fout.write(f'{nellipsoids} ellipsoids\n')
    fout.write('\n')

    # Box dimensions ----------------------------------------------------------
    fout.write('{} {} xlo xhi\n'.format(xlo, xhi))
    fout.write('{} {} ylo yhi\n'.format(ylo, yhi))
    fout.write('{} {} zlo zhi\n'.format(zlo, zhi))
    fout.write('\n')

    # Atoms section -----------------------------------------------------------
    fout.write('Atoms # hybrid\n\n')

    jadd = 1
    for j in range(natoms0):
        fout.write('{:.0f} {:.0f} {:.6f} {:.6f} {:.6f}'\
                   ' {:.0f} {:.4f} {:.0f} {:.0f} 0 0 0'\
                   '\n'.format(j+jadd,types0,*tuple(atoms0[j,:]),
                   ellipsoidflag0,density0,molecule_ID0,charge0))

    fout.write('\n')

    # Ellipsoids section ------------------------------------------------------
    fout.write('Ellipsoids\n\n');

    jadd = 1
    for j in range(natoms0):
            fout.write('{:.0f} {:.2f} {:.2f} {:.2f} {:.6f} '\
                       '{:.6f} {:.6f} {:.6f}\n'.format(j+jadd,
                       *tuple(shape0),*tuple(quat0[j,:])))

    fout.write('\n')

print('init.data generated.')

# Exit ========================================================================
sys.exit()
