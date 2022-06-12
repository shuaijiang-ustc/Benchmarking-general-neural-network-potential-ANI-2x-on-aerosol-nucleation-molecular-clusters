import os
import torch
import torchani
import ase.io
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt

def _sort(list,a,b):

    list.sort(key = lambda x:int(x.split(a)[0].split(b)[-1]))
    return list

Hatree_to_kcal_per_mol=627.5094
eV_to_Hatree=1/27.2114

method_dirs=['ANI_opted_structures_energies', 'PM7_opted_structures_energies']
cur_dir=os.getcwd()

E_DFT=[]
E_ANI=[]
E_PM7=[]

target_path=os.path.join(cur_dir, 'DFT_opted_structures_energies')
os.chdir(target_path)
all_files=os.listdir(target_path)
all_xyzs=[]
for file_name in all_files:
	if file_name[-4:]=='.xyz':
		all_xyzs.append(file_name)
all_xyzs_sorted=_sort(all_xyzs, '.', '_')
for dft_xyz in all_xyzs_sorted:
	with open(dft_xyz,'r') as f:
		f.readline()
		energy_line=f.readline()
		E=float(energy_line.split()[0])
		E_DFT.append(E)		
for method_dir in method_dirs:
	target_path=os.path.join(cur_dir, method_dir)
	os.chdir(target_path)
	for dft_xyz in all_xyzs_sorted:
		if 'ANI' in method_dir:
			xyz='ANI_optimized_Geom_'+dft_xyz.split('.')[0].split('_')[-1]+'.xyz'
			with open(xyz,'r') as f:
				f.readline()
				energy_line=f.readline()
				E=float(energy_line.split('free_energy=')[0].split('energy=')[-1])*eV_to_Hatree
				E_ANI.append(E)
		else:
			xyz='Geom_'+dft_xyz.split('.')[0].split('_')[-1]+'.xyz'
			#print(xyz)		
			with open(xyz,'r') as f:
				f.readline()
				energy_line=f.readline()
				E=float(energy_line.split()[0])
				E_PM7.append(E*eV_to_Hatree)
xyz_atoms_objs = ase.io.read('%s'%xyz, index=':')
composition=xyz_atoms_objs[0].get_atomic_numbers().tolist()
os.chdir(cur_dir)
E_DFT=np.array(E_DFT) # in Hatree
E_ANI=np.array(E_ANI) # in Hatree
E_PM7=np.array(E_PM7) # in Hatree

# obtain monomers E
monomer_set=[]
S=0
N=0
C=0
O=0
for element in composition:
	if element==16:
		S+=1
	elif element==6:
		C+=1
	elif element==7:
		N+=1
	elif element==8:
		O+=1
if S==2 and O==8:
	monomer_set.append('SA')
	monomer_set.append('SA')
elif S==1 and O==4 and C==2 and N==1:
	monomer_set.append('SA')
	monomer_set.append('DMA')
elif S==1 and O==5:
	monomer_set.append('SA')
	monomer_set.append('W')	
elif S==1 and O==4 and N==1 and C==0:
	monomer_set.append('SA')
	monomer_set.append('Am')
elif S==1 and O==8 and C==4:
	monomer_set.append('SA')
	monomer_set.append('SUA')
i=0	
for monomer in monomer_set:
	i+=1
	if i==1:
		f=open('../../monomers_optimized/DFT_optimized_%s.xyz'%monomer,'r')
		f.readline()
		energy_line=f.readline()
		E_monomer_A_DFT=float(energy_line.split()[0])
		f.close()
		
		f=open('../../monomers_optimized/ANI_optimized_%s.xyz'%monomer,'r')
		f.readline()
		energy_line=f.readline()
		E_monomer_A_ANI=float(energy_line.split('free_energy=')[0].split('energy=')[-1])*eV_to_Hatree
		f.close()
		
		f=open('../../monomers_optimized/%s.arc'%monomer,'r')
		for line in f.readlines():
			if 'TOTAL ENERGY' in line:
				E_monomer_A_pm7=float(line.split()[-2])*eV_to_Hatree
		f.close()
	elif i==2:
		f=open('../../monomers_optimized/DFT_optimized_%s.xyz'%monomer,'r')
		f.readline()
		energy_line=f.readline()
		E_monomer_B_DFT=float(energy_line.split()[0])
		f.close()
		
		f=open('../../monomers_optimized/ANI_optimized_%s.xyz'%monomer,'r')
		f.readline()
		energy_line=f.readline()
		E_monomer_B_ANI=float(energy_line.split('free_energy=')[0].split('energy=')[-1])*eV_to_Hatree
		f.close()
		
		f=open('../../monomers_optimized/%s.arc'%monomer,'r')
		for line in f.readlines():
			if 'TOTAL ENERGY' in line:
				E_monomer_B_pm7=float(line.split()[-2])*eV_to_Hatree
		f.close()
BE_DFT=[]
for DFT in list(E_DFT):
	BE=DFT-E_monomer_A_DFT-E_monomer_B_DFT
	BE_DFT.append(BE)
BE_DFT=np.array(BE_DFT)

E_ANI=list(E_ANI)
BE_ANI=[]
for E in E_ANI:
	BE_ANI.append(E-E_monomer_A_ANI-E_monomer_B_ANI)
BE_ANI=np.array(BE_ANI)

BE_PM7=[]
for PM7 in list(E_PM7):
	BE=PM7-E_monomer_A_pm7-E_monomer_B_pm7
	BE_PM7.append(BE)
BE_PM7=np.array(BE_PM7)

BE_DFT=BE_DFT*Hatree_to_kcal_per_mol
BE_ANI=BE_ANI*Hatree_to_kcal_per_mol
BE_PM7=BE_PM7*Hatree_to_kcal_per_mol

f=open('ANI_DFT_PM7_BE.txt','w')
f.write('BE_DFT (kcal/mol)'+'\t'+'BE_ANI (kcal/mol)'+'\t'+'BE_PM7 (kcal/mol)'+'\n')
for i in range(len(BE_ANI)):
	BE_DFT_kcal_mol=format(BE_DFT[i],'.2f')
	BE_ANI_kcal_mol=format(BE_ANI[i],'.2f')
	BE_PM7_kcal_mol=format(BE_PM7[i],'.2f')
	f.write(BE_DFT_kcal_mol+'\t'+BE_ANI_kcal_mol+'\t'+BE_PM7_kcal_mol+'\n')
f.close()

