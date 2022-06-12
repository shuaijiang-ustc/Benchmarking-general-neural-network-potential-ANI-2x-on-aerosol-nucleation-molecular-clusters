
import ase.io
import numpy as np
import ase
import matplotlib.pylab as plt

xyz_file_name='./SA_DMA/bowl_0.05/right_CNs_selected.xyz'

xyz_atoms_objs = ase.io.read('%s'%xyz_file_name, index=':')
proton_trans_set=[]
noPT=[]
noPT_noHB=[]
noPT_oneHB=[]
noPT_twoHB=[]
for i in range(len(xyz_atoms_objs)):
	N_H_bond_numb=0
	N_HB=0
	species=xyz_atoms_objs[0].get_atomic_numbers()
	index = np.argwhere(species == 7)
	for index_other, element in np.ndenumerate(species):
		if species[index_other]==1: # H
			N_H_distance=xyz_atoms_objs[i].get_distance(index, index_other)
			if N_H_distance< 1.27788888405839:
				N_H_bond_numb+=1
				for index_other1, element1 in np.ndenumerate(species):
					if species[index_other1]==4: # O
						H_O_distance=xyz_atoms_objs[i].get_distance(index_other, index_other1)
						if H_O_distance > 1.23356069976309 and H_O_distance < 2.0: # N-H...O
							N_HB+=1
			elif N_H_distance > 1.27788888405839 and N_H_distance < 2.0: # N...H
				N_HB+=1
	if N_H_bond_numb==2: # two N-H bonds
		proton_trans_set.append(i)
	else:
		noPT.append(i)
		if N_HB==0:
			noPT_noHB.append(i)
		elif N_HB==1:
			noPT_oneHB.append(i)
		elif N_HB==2:
			noPT_twoHB.append(i)

with open('./SA_DMA/bowl_0.05/ANI_DFT_RE.txt', 'r') as f:
	f.readline()
	lines=f.readlines()
	
delta_RE_PT=[]
delta_RE_noPT=[]
delta_RE_noPT_noHB=[]
delta_RE_noPT_oneHB=[]
delta_RE_noPT_twoHB=[]
for i, line in enumerate(lines):
	RE_DFT=float(line.split()[0])
	RE_ANI=float(line.split()[1])
	delta_RE=RE_ANI-RE_DFT
	if i in proton_trans_set:
		delta_RE_PT.append(delta_RE)
	if i in noPT:
		delta_RE_noPT.append(delta_RE)
	if i in noPT_noHB:
		delta_RE_noPT_noHB.append(delta_RE)
	if i in noPT_oneHB:
		delta_RE_noPT_oneHB.append(delta_RE)
	if i in noPT_twoHB:
		delta_RE_noPT_twoHB.append(delta_RE)

x=np.array(['PT',  'no_PT_no_HB', 'no_PT_one_HB', 'no_PT_two_HBs'])
y=np.array([np.mean(delta_RE_PT),  np.mean(delta_RE_noPT_noHB), np.mean(delta_RE_noPT_oneHB), np.mean(delta_RE_noPT_twoHB)])
y_min=np.array([min(delta_RE_PT),  min(delta_RE_noPT_noHB), min(delta_RE_noPT_oneHB), min(delta_RE_noPT_twoHB)])
y_max=np.array([max(delta_RE_PT),  max(delta_RE_noPT_noHB), max(delta_RE_noPT_oneHB), max(delta_RE_noPT_twoHB)])
yerr = np.zeros([2,len(y)])
yerr[0,:] = y - y_min
yerr[1,:] = y_max - y
plt.ylabel(r"RE_ANI-2x - RE_DFT (kcal/mol)")
plt.errorbar(x,y,yerr=yerr[:,:],ecolor='k',elinewidth=0.5,marker='s',mfc='orange',\
	mec='k',mew=1,ms=10,alpha=1,capsize=5,capthick=3,linestyle="none")
plt.savefig('fig5_RE_PT_HB_correlation.jpeg', dpi=500)
#plt.show()

		
		
			