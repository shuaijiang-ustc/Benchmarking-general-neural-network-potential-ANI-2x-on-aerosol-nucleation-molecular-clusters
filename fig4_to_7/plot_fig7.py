import os
from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.metrics import mean_squared_error
from matplotlib.pyplot import MultipleLocator
import matplotlib.gridspec as gridspec

KJPERHARTREE = 2625.499638
JOULEPERHARTREE = KJPERHARTREE*1000.0
Hatree_to_ev=27.2114

color_other=['#4456c7', '#f27d51', '#fdc765','#b5a6a7']
dpi_numb=600
linewidth_value=2.0
markersize_value=4
fontsize=8
x_major_locator=MultipleLocator(10)
y_major_locator=MultipleLocator(10)
a=11
b=4
tick_label_size=16
font_legend = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 11,
}
font_title = {
'weight' : 'normal',
'size'   : 18,
}


cur_dir=os.getcwd()
dir_list=['SA_W\\bowl_0.05', 'SA_Am\\bowl_0.05', 'SA_DMA\\bowl_0.05', 'SA_SA\\bowl_0.05', 'SA_SUA\\bowl_0.05']
system_list=['(SA)$_{1}$(W)$_{1}$', '(SA)$_{1}$(Am)$_{1}$','(SA)$_{1}$(DMA)$_{1}$','(SA)$_{2}$','(SA)$_{1}$(SUA)$_{1}$' ]

def square_with_sign(input):
	if input<0.0:
		output=input**2*(-1)
	else:
		output=input**2
	return output
def sqrt_with_sign(input):
	if input<0.0:
		output=sqrt(input*(-1))*(-1)
	else:
		output=sqrt(input)
	return output

def get_force(force_file_name, force_type):	
	if force_type=='DFT' or 'PM7' or 'ANI': # force unit is J/mol, should be converted to eV/A
		f=open(force_file_name,'r')
		force=[]
		for line in f.readlines():
			if len(line.split())==3:
				force_x_single=float(line.split()[0])/JOULEPERHARTREE*Hatree_to_ev
				force_y_single=float(line.split()[1])/JOULEPERHARTREE*Hatree_to_ev
				force_z_single=float(line.split()[2])/JOULEPERHARTREE*Hatree_to_ev
				force_single=sqrt_with_sign(square_with_sign(force_x_single)+square_with_sign(force_y_single)+square_with_sign(force_z_single))
				force.append(force_single)
		f.close()
	return force

i=0
ax_set=[]
fig=plt.figure(num=1, figsize=(8.5, 5.5), dpi=600, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(2, 6)
for dirt in dir_list:
	i+=1
	target_path=os.path.join(cur_dir, dirt)
	target_file=os.path.join(target_path, 'Force_J_per_mol_Nuclear_charge.txt')
	DFT_force=get_force(target_file, 'DFT')
	target_file=os.path.join(target_path, 'Force_ANI_J_per_mol_ANI.txt')
	ANI_force=get_force(target_file, 'ANI')
	target_file=os.path.join(target_path, 'Force_J_per_mol_PM7.txt')
	PM7_force=get_force(target_file, 'PM7')
	
	rmse_ANI = format(sqrt(mean_squared_error(DFT_force, ANI_force)),'.2f')
	rmse_PM7 = format(sqrt(mean_squared_error(DFT_force, PM7_force)),'.2f')
	
	if i==1:
		ax = plt.subplot(gs[0, :2])
	elif i==2:
		ax = plt.subplot(gs[0, 2:4])
	elif i==3:
		ax = plt.subplot(gs[0, 4:6])
	elif i==4:
		ax = plt.subplot(gs[1, 1:3])
	elif i==5:
		ax = plt.subplot(gs[1, 3:5])
	ax_set.append(ax)
	plt.scatter(DFT_force,  ANI_force, marker="o", color=color_other[0], s=markersize_value)
	plt.scatter(DFT_force,  PM7_force, marker="o", color=color_other[1], s=markersize_value)
	ax.xaxis.set_major_locator(x_major_locator)
	ax.yaxis.set_major_locator(y_major_locator)
	axis_low=min(min(DFT_force), min(ANI_force), min(PM7_force))
	axis_high=max(max(DFT_force), max(ANI_force), max(PM7_force))
	axis_low_about=(axis_low//10)*10
	axis_high_about=(axis_high//10+1)*10
	x=np.arange(axis_low_about,axis_high_about,(axis_high_about-axis_low_about)/100.0)
	y=x
	plt.plot(x,y, color=color_other[3], linewidth=linewidth_value)
	plt.title(system_list[i-1])
	plt.xlim(axis_low_about, axis_high_about)
	plt.ylim(axis_low_about, axis_high_about)
	ANGSTROM= "Å"
	labels = ax.get_xticklabels() + ax.get_yticklabels()
	plt.ylabel('ANI-2x/PM7 Force (eV/%s)'%ANGSTROM)
	plt.xlabel('DFT Force (eV/%s)'%ANGSTROM)
	ax.text(0.015,0.93, 'RMSE_ANI-2x=%s (eV/%s)'%(rmse_ANI,ANGSTROM), transform=ax.transAxes,fontsize=fontsize)
	ax.text(0.015,0.85, 'RMSE_PM7=%s (eV/%s)'%(rmse_PM7,ANGSTROM), transform=ax.transAxes,fontsize=fontsize)
	
fig.legend(ax_set, labels=["ANI-2x", "PM7"],
           loc="upper right")
plt.tight_layout()	
plt.savefig('./fig7_Force_on_MetaMD.jpeg',dpi=dpi_numb)	                                                                                                                                     