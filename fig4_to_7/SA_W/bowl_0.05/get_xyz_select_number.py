import ase.io
import math

select_numb=500


xyz_atoms_objs=ase.io.read('right_CNs.xyz',index=':')

frame_numb=len(xyz_atoms_objs)
print("total xyz number before selection: %d"%frame_numb)
every_numb=math.floor(frame_numb/select_numb)

i_selected=[]
for i in range(frame_numb):
		if (i+1)%every_numb==0:
			i_selected.append(i)
		
ase.io.write('./right_CNs_selected.xyz', [xyz_atoms_objs[i] for i in i_selected[:select_numb]])

xyz_atoms_objs1=ase.io.read('right_CNs_selected.xyz',index=':')
print("total xyz number after selection: %d"%len(xyz_atoms_objs1))
