# -*- coding: utf-8 -*-
"""
Structure minimization and constant temperature MD using ASE interface
======================================================================

This example is modified from the official `home page` and
`Constant temperature MD`_ to use the ASE interface of TorchANI as energy
calculator.

.. _home page:
    https://wiki.fysik.dtu.dk/ase/
.. _Constant temperature MD:
    https://wiki.fysik.dtu.dk/ase/tutorials/md/md.html#constant-temperature-md
"""


###############################################################################
# To begin with, let's first import the modules we will use:
from ase.optimize import BFGS
from ase import units
import torchani
import ase.io
import os

calculator = torchani.models.ANI2x().ase()
for file_name in os.listdir('../PM7_opted_structures_energies/'):
	if file_name[-4:]=='.xyz':
		xyz_atoms_objs=ase.io.read('../PM7_opted_structures_energies/%s'%file_name,index=':')
		xyz_atoms_objs[0].set_calculator(calculator)
		print("Begin minimizing...")
		opt = BFGS(xyz_atoms_objs[0])
		opt.run(fmax=0.001)
		output_xyz_name='ANI_optimized_'+file_name
		ase.io.write('./%s'%output_xyz_name, xyz_atoms_objs[0])
