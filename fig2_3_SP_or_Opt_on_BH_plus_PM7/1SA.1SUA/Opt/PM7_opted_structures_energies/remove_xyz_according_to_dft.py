

# obtain the isomers names in DFT dir

import os


dft_dir='../DFT_opted_structures_energies'
isomers_index=[]
for file_name in os.listdir(dft_dir):
	if 'Geom' in file_name:
		isomers_index.append(file_name[:-4].split('_')[-1])

for file_name in os.listdir(os.getcwd()):
	if 'Geom' in file_name:
		isomer_index=file_name[:-4].split('_')[-1]
		if isomer_index in isomers_index:
			pass
		else:
			os.system('rm %s'%file_name)


#print(isomers_index)
#print(len(isomers_index))
		