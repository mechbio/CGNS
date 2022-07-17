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

# date         :11-Jul-22
# version      :0.9.0
# usage        :python3 ne.py 2.5
# py_version   :3.8.10

# Initialize ==================================================================

import sys
import numpy as np
import random

scl = float(sys.argv[1]) # scale of simulation box
L = 120*scl # Length of sides of simulation box in x and y (Lx = Ly = L)
D0 = 100*scl # Diameter of spherical_monolayer
sd = 221 # seed for random generator

mNM = 1.237
mNPC = 60.00
mNLA = 0.070
mNLB = 0.068

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

    for i in range(1,int(np.pi/2.0/np.arcsin(d/2.0/R1))):
        ti = (i-1)*np.pi/int(np.pi/2.0/np.arcsin(d/2.0/R1))
        for j in range(0,int(np.pi*np.sin(ti)/np.arcsin(d/2.0/R1))):
            fj = (j-1)*2.0*np.pi/int(np.pi*np.sin(ti)/np.arcsin(d/2.0/R1))
            x0 = R1*np.sin(ti)*np.cos(fj)
            y0 = R1*np.sin(ti)*np.sin(fj)
            z0 = R1*np.cos(ti)
            xyz = np.vstack([xyz, [x0, y0, z0]])

    xyz = xyz + np.random.rand(np.shape(xyz)[0],3)/100

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

def distrib(Nb,NBf,NAf,ft0,ft1):
    """Filament distribution.
    """

    Bf = []
    Af = []
    for i in range(len(Nb)):
        Bf = np.concatenate((Bf, Nb[i]*np.ones(NBf[i])))
        Af = np.concatenate((Af, Nb[i]*np.ones(NAf[i])))

    nfil = np.sum(NBf)+np.sum(NAf)
    orderedtypes = np.concatenate((ft0*np.ones(np.sum(NAf)),
                                   ft1*np.ones(np.sum(NBf))))
    orderedNbs = np.concatenate((Af, Bf))
    p = np.random.RandomState(seed=sd).permutation(nfil)
    filtype = orderedtypes[p].astype(int)
    filNb = orderedNbs[p].astype(int)

    return filtype, filNb

def spherical_filaments(D,d,filtype,filNb,natoms,nbonds,nangles,nmols):
    """Filaments spaced by d grown radially inwards in a sphere of diameter D.
    """
    atoms = np.zeros(9)
    bonds = np.zeros(4)
    angles = np.zeros(5)

    R1 = D/2
    # xyz = np.asarray([[0.0, 0.0, R1],
    #                     [0.0, 0.0, -R1]])

    nfil = len(filtype)
    atom_id = natoms
    bond_id = nbonds
    angle_id = nangles
    fil = 0
    for i in range(1,int(np.pi/2.0/np.arcsin(d/2.0/R1))):
        ti = (i-1)*np.pi/int(np.pi/2.0/np.arcsin(d/2.0/R1))
        for j in range(0,int(np.pi*np.sin(ti)/np.arcsin(d/2.0/R1))):
            fj = (j-1)*2.0*np.pi/int(np.pi*np.sin(ti)/np.arcsin(d/2.0/R1))
            if fil < nfil:
                fil += 1
                seq=0
                Nb = filNb[fil-1]
                for k in range(Nb):
                    seq += 1
                    atom_id += 1

                    atomtype = filtype[fil-1]

                    if filtype[fil-1] == 3:
                        btype = 1
                        atype = 1
                        atommass = mNLA
                        ri = -2.0
                    if filtype[fil-1] == 4:
                        btype = 1
                        atype = 1
                        atommass = mNLB
                        ri = -2.0

                    x0 = (R1+ri*k)*np.sin(ti)*np.cos(fj) + np.random.rand()/100
                    y0 = (R1+ri*k)*np.sin(ti)*np.sin(fj) + np.random.rand()/100
                    z0 = (R1+ri*k)*np.cos(ti) + np.random.rand()/100
                    atoms0 = [atom_id, atomtype, x0, y0,
                              z0, 0, atommass,
                              nmols+fil, k]
                    atoms = np.vstack([atoms, atoms0])

                    if seq<Nb:
                        bond_id += 1
                        bonds0 = [bond_id, btype, atom_id, atom_id+1]
                        bonds = np.vstack([bonds, bonds0])

                    if seq>1 and seq<Nb:
                        angle_id += 1
                        angles0 = [angle_id, atype,
                                   atom_id-1, atom_id, atom_id+1]
                        angles = np.vstack([angles, angles0])

    atoms = atoms[1:]
    bonds = bonds[1:]
    angles = angles[1:]

    return atoms, bonds, angles, fil

# Molecules ===================================================================

# NM --------------------------------------------------------------------------
d0 = 2.0
d1 = d0*1.08 # NM collapses a bit
d1 = d0*0.97

natoms0, quat0, atoms0 = spherical_monolayer(D0,d1)
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

nNPC = int(6*sclsq*np.pi)
listNPC = np.ones(nNPC)
nNM = int(natoms0-nNPC)
listNM = np.ones(nNM)
orderedtypes = np.concatenate((2*listNPC, 1*listNM))
orderedmasses = np.concatenate((mNPC*listNPC, mNM*listNM))
p = np.random.RandomState(seed=sd).permutation(natoms0)
types0 = orderedtypes[p]
mass = orderedmasses[p]
shape0 = [d0, d0, d0]
volm = (4.0/3.0*(np.pi*shape0[0]*shape0[1]*shape0[2])/8.0)
density0 = mass/volm
ellipsoidflag0 = 1
molecule_ID0 = 0
charge0 = 0

natoms += natoms0
nellipsoids += natoms0
natomTs += 2

print(f'{natoms} NM particles generated.')

# NL --------------------------------------------------------------------------
d1 = 8.14
h = -5.0
Nb  = [9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 43]
NBf = np.asarray([3, 12, 16, 11, 11, 7, 4, 2, 1, 0, 1, 0])*sclsq*np.pi
NAf = np.asarray([7, 24, 33, 21, 23, 15, 8, 3, 3, 1, 0, 1])*sclsq*np.pi
NAf = np.rint(NAf).astype(int)
NBf = np.rint(NBf).astype(int)

(filtype, filNb) = distrib(Nb,NBf,NAf,3,4)
(atoms1, bonds1, angles1, fil) = spherical_filaments(D0-2.0,d1,filtype,filNb,
                                 natoms,nbonds,nangles,nmols)

natoms1 = np.shape(atoms1)[0]
nbonds1 = np.shape(bonds1)[0]
nangles1 = np.shape(angles1)[0]

natoms += natoms1
nbonds += nbonds1
nangles += nangles1
natomTs += 2 # NLA, NLB
nbondTs += 3 # NL fil, CL, IL
nangleTs += 1 # NL fil
nmols += fil

print(f'{fil} of {len(filtype)} NL filaments generated.')

fracf = np.around(100*(NAf+NBf)/np.sum(NAf+NBf),2)
print(f'NL filament distribution: {fracf}.')

# Write =======================================================================

with open('ne.data','w') as fout:
    # Title comment -----------------------------------------------------------
    fout.write('Initial configuration of NE\n\n')

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
                   '\n'.format(j+jadd,types0[j],*tuple(atoms0[j,:]),
                   ellipsoidflag0,density0[j],molecule_ID0,charge0))

    jadd += natoms0
    for j in range(natoms1):
        fout.write('{:.0f} {:.0f} {:.6f} {:.6f} {:.6f}'\
                   ' {:.0f} {:.4f} {:.0f} {:.0f} 0 0 0'\
                   '\n'.format(*tuple(atoms1[j,:])))

    fout.write('\n')

    # Bonds section -----------------------------------------------------------
    fout.write('Bonds\n\n')

    jadd = 1
    for j in range(nbonds1):
        fout.write('{:.0f} {:.0f} {:.0f} {:.0f}'\
                   '\n'.format(*tuple(bonds1[j,:])))

    fout.write('\n')

    # Angles section ----------------------------------------------------------
    fout.write('Angles\n\n')

    jadd = 1
    for j in range(nangles1):
        fout.write('{:.0f} {:.0f} {:.0f} {:.0f} {:.0f}'\
                   '\n'.format(*tuple(angles1[j,:])))


    fout.write('\n')

    # Ellipsoids section ------------------------------------------------------
    fout.write('Ellipsoids\n\n');

    jadd = 1
    for j in range(natoms0):
            fout.write('{:.0f} {:.2f} {:.2f} {:.2f} {:.6f} '\
                       '{:.6f} {:.6f} {:.6f}\n'.format(j+jadd,
                       *tuple(shape0),*tuple(quat0[j,:])))

    fout.write('\n')

print('ne.data generated.')

# Exit ========================================================================
sys.exit()
