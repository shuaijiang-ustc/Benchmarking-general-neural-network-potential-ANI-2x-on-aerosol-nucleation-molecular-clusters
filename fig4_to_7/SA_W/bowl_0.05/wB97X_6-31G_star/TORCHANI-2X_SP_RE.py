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
xyz_file_name='right_CNs_selected.xyz'

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model = torchani.models.ANI2x(periodic_table_index=True).to(device)
xyz_atoms_objs = ase.io.read('%s'%xyz_file_name, index=':')
for i in range(len(xyz_atoms_objs)):
	coordinates = torch.from_numpy(xyz_atoms_objs[i].get_positions()).type(torch.float)
	coordinates=coordinates.reshape(1,-1,3)
	coordinates.requires_grad_()  
	species=torch.from_numpy(xyz_atoms_objs[i].get_atomic_numbers()).type(torch.long)
	composition=species.tolist()
	species=species.reshape(1,-1)
	energy = model((species, coordinates)).energies
	derivative = torch.autograd.grad(energy.sum(), coordinates)[0]
	force = -derivative
	if i==0:
		tmp_E=energy
	else:
		tmp_E=torch.cat([tmp_E, energy], dim=0)
	print('Index: %d, Energy:'%i, energy.item())

E_ANI=tmp_E.detach().numpy()
f=open('Energy_Hatree.txt','r')
E_DFT=[]
for line in f.readlines():
	E_DFT.append(float(line.split()[0]))
f.close()
E_DFT=np.array(E_DFT)

# caculate relative energy 
RE_DFT=[]
RE_ANI=[]
E_DFT_min=min(E_DFT)
E_ANI_min=min(E_ANI)
for i in range(len(E_DFT)):
	RE_DFT.append(E_DFT[i]-E_DFT_min)
	RE_ANI.append(E_ANI[i]-E_ANI_min)
RE_DFT=np.array(RE_DFT)
RE_ANI=np.array(RE_ANI)
RE_DFT=RE_DFT*Hatree_to_kcal_per_mol
RE_ANI=RE_ANI*Hatree_to_kcal_per_mol
f=open('ANI_DFT_RE.txt','w')
f.write('RE_DFT (kcal/mol)'+'\t'+'RE_ANI (kcal/mol)'+'\n')
for i in range(len(RE_ANI)):
	RE_DFT_kcal_mol=format(RE_DFT[i],'.2f')
	RE_ANI_kcal_mol=format(RE_ANI[i],'.2f')
	f.write(RE_DFT_kcal_mol+'\t'+RE_ANI_kcal_mol+'\n')
f.close()