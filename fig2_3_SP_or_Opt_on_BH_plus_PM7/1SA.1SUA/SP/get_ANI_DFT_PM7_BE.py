import torch
import torchani
import ase.io
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt

color_other=['#4456c7', '#f27d51', '#fdc765','#b5a6a7']
dpi_numb=600
linewidth_value=2.0
markersize_value=4
fontsize=16
#x_major_locator=MultipleLocator(10)
#y_major_locator=MultipleLocator(10)
a=11
b=4
tick_label_size=16
font_legend = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 11,
}
font_title = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 18,
}


Hatree_to_kcal_per_mol=627.5094
xyz_file_name='../BH/combined_xyz.xyz'
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = torchani.models.ANI2x(periodic_table_index=True).to(device)
xyz_atoms_objs = ase.io.read('%s'%xyz_file_name, index=':')
coordinates = torch.from_numpy(xyz_atoms_objs[0].get_positions()).type(torch.float)
coordinates=coordinates.reshape(1,-1,3)
coordinates.requires_grad_()  
species=torch.from_numpy(xyz_atoms_objs[0].get_atomic_numbers()).type(torch.long)
composition=species.tolist()
species=species.reshape(1,-1)

f=open('./DFT/Energy_Hatree.txt','r')
E_DFT=[]
for line in f.readlines():
	E_DFT.append(float(line.split()[0]))
f.close()
E_DFT=np.array(E_DFT)

f=open('./ANI/Energy_Hatree_ANI.txt','r')
E_ANI=[]
for line in f.readlines():
	E_ANI.append(float(line.split()[0]))
f.close()
E_ANI=np.array(E_ANI)

f=open('./PM7/Energy_Hatree_PM7.txt','r')
E_PM7=[]
for line in f.readlines():
	E_PM7.append(float(line.split()[0]))
f.close()
E_PM7=np.array(E_PM7)

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
		f=open('../../monomers/%s.xyz'%monomer,'r')
		f.readline()
		energy_line=f.readline()
		E_monomer_A=float(energy_line.split()[0])
		E_monomer_A_PM7=float(energy_line.split()[2])
		f.close()
		xyz_monomer_A=ase.io.read('../../monomers/%s.xyz'%monomer, index=':')
		coordinates=torch.from_numpy(xyz_monomer_A[0].get_positions()).type(torch.float)
		coordinates=coordinates.reshape(1,-1,3)
		coordinates.requires_grad_()
		species=torch.from_numpy(xyz_monomer_A[0].get_atomic_numbers()).type(torch.long)
		species=species.reshape(1,-1)
		energy_A = model((species, coordinates)).energies 
	elif i==2:
		f=open('../../monomers/%s.xyz'%monomer,'r')
		f.readline()
		energy_line=f.readline()
		E_monomer_B=float(energy_line.split()[0])
		E_monomer_B_PM7=float(energy_line.split()[2])
		f.close()
		xyz_monomer_B=ase.io.read('../../monomers/%s.xyz'%monomer, index=':')
		coordinates=torch.from_numpy(xyz_monomer_B[0].get_positions()).type(torch.float)
		coordinates=coordinates.reshape(1,-1,3)
		coordinates.requires_grad_()
		species=torch.from_numpy(xyz_monomer_B[0].get_atomic_numbers()).type(torch.long)
		species=species.reshape(1,-1)
		energy_B = model((species, coordinates)).energies
BE_DFT=[]
for DFT in list(E_DFT):
	BE=DFT-E_monomer_A-E_monomer_B
	BE_DFT.append(BE)
BE_DFT=np.array(BE_DFT)

BE_PM7=[]
for PM7 in list(E_PM7):
	BE=PM7-E_monomer_A_PM7-E_monomer_B_PM7
	BE_PM7.append(BE)
BE_PM7=np.array(BE_PM7)

E_ANI=list(E_ANI)
BE_ANI=[]
for E in E_ANI:
	BE_ANI.append(E-energy_A.item()-energy_B.item())
BE_ANI=np.array(BE_ANI)

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