import ase.io
import os
import numpy as np

KJPERHARTREE = 2625.499638
JOULEPERHARTREE = KJPERHARTREE*1000.0
BOHRPERA = 1.889725989

def G16_PM7(crds, atoms, basis_ = '6-31++G(d,p)', xc_='wB97XD', jobtype_='force', threads=16):
	istring = ""
	for j in range(len(atoms)):
		istring=istring+str(atoms[j])+' '+str(crds[j,0])+' '+str(crds[j,1])+' '+str(crds[j,2])+'\n'
	line2='''%%nproc=%d        
%%mem=2000MB        
#p %s %s Force 
                     
aaaaa                
                     
0 1
%s
'''%(threads,xc_,basis_,istring)
	
	current_dir=os.getcwd()
	
	ff=open("init.gjf","w")
	ff.write(line2)
	ff.close()
	os.system("g16 %s %s"%("init.gjf","init.log"))
	
	f=open("init.log")
	lines=f.readlines()
	f.close()
	Energy=0                              
	Forces = np.zeros((len(crds),3))
	for i, line in enumerate(lines):
		if line.count('SCF Done:')>0:
			Energy = float(line.split()[4])
		if line.count('Forces (Hartrees/Bohr)') > 0:
			k = 0
			for j in range(1, len(crds)+1):
				Forces[j-1,:] = float(lines[i+k+3].split()[2])*(-1), float(lines[i+k+3].split()[3])*(-1), float(lines[i+k+3].split()[4])*(-1)
				k += 1
			break
	f1=open("Energy_Hatree_PM7.txt","a+")
	f1.write(str(Energy))
	f1.write('\n')
	f1.close()	
	
	Forces_unit=-Forces*JOULEPERHARTREE*BOHRPERA# Hatree/Bohr => J/mol
	with open('Force_J_per_mol_PM7.txt', 'a') as f:
		f.write(str(len(atoms))+'\n')
		for item in atoms:
			f.write(str(item)+' ')
		f.write('\n')
	f.close()
	with open('Force_J_per_mol_PM7.txt', 'ab') as f:
		np.savetxt(f, Forces_unit)
	f.close() 
	
	return Energy, -Forces*JOULEPERHARTREE*BOHRPERA# Hatree/Bohr => J/mol
	

xyz_atoms_objs=ase.io.read('right_CNs_selected.xyz',index=':')
for i in range(len(xyz_atoms_objs)):  
	G16_PM7(xyz_atoms_objs[i].get_positions(), xyz_atoms_objs[i].get_atomic_numbers(), basis_ = '', xc_='PM7', jobtype_='force', threads=8)