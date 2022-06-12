from ase.optimize import BFGS
from ase import units
import torchani
import ase.io
import os

calculator = torchani.models.ANI2x().ase()
for file_name in os.listdir('.'):
	if file_name[-4:]=='.xyz':
		xyz_atoms_objs=ase.io.read(file_name,index=':')
		xyz_atoms_objs[0].set_calculator(calculator)
		print("Begin minimizing...")
		opt = BFGS(xyz_atoms_objs[0])
		opt.run(fmax=0.001)
		output_xyz_name='ANI_optimized_'+file_name
		ase.io.write('./%s'%output_xyz_name, xyz_atoms_objs[0])
