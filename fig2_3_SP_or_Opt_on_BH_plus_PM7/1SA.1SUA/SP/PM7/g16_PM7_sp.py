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
#p %s %s sp 
                     
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
	#Dipole=[]
	for i, line in enumerate(lines):
		if line.count('SCF Done:')>0:
			Energy = float(line.split()[4])
	f1=open("Energy_Hatree_PM7.txt","a+")
	f1.write(str(Energy))
	f1.write('\n')
	f1.close()	

xyz_atoms_objs=ase.io.read('../../BH/combined_xyz.xyz',index=':')
for i in range(len(xyz_atoms_objs)):  
	G16_PM7(xyz_atoms_objs[i].get_positions(), xyz_atoms_objs[i].get_atomic_numbers(), basis_ = '', xc_='PM7', jobtype_='force', threads=8)