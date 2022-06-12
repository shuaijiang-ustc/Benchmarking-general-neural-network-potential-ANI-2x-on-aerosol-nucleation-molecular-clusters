import os
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
from matplotlib.pyplot import MultipleLocator
import matplotlib.gridspec as gridspec

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
font_title = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 18,
}
cur_dir=os.getcwd()
dir_list=['1SA.1H2O', '1SA.1NH3', '1SA.1DMA', '2SA', '1SA.1SUA']
system_list=['(SA)$_{1}$(W)$_{1}$', '(SA)$_{1}$(Am)$_{1}$','(SA)$_{1}$(DMA)$_{1}$','(SA)$_{2}$','(SA)$_{1}$(SUA)$_{1}$']


i=0
ax_set=[]
fig1=plt.figure(num=1, figsize=(8.5, 5.5), dpi=600, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(2, 6)
for dirt in dir_list:
	i+=1
	target_path=os.path.join(cur_dir, dirt)
	target_path=os.path.join(target_path, 'SP')
	target_file=os.path.join(target_path, 'ANI_DFT_PM7_BE.txt')
	f=open(target_file,'r')
	f.readline()
	BE_DFT=[]
	BE_ANI=[]
	BE_PM7=[]
	for line in f.readlines():
		BE_DFT.append(float(line.split()[0]))
		BE_ANI.append(float(line.split()[1]))
		BE_PM7.append(float(line.split()[2]))
	f.close()
	RMSE_BE = metrics.mean_squared_error(BE_ANI, BE_DFT)**0.5
	RMSE_BE_PM7 = metrics.mean_squared_error(BE_PM7, BE_DFT)**0.5
	fig1=plt.figure(num=1, figsize=(8, 6), dpi=600, facecolor='w', edgecolor='k')
	if i==1:
		ax1 = plt.subplot(gs[0, :2])
	elif i==2:
		ax1 = plt.subplot(gs[0, 2:4])
	elif i==3:
		ax1 = plt.subplot(gs[0, 4:6])
	elif i==4:
		ax1 = plt.subplot(gs[1, 1:3])
	elif i==5:
		ax1 = plt.subplot(gs[1, 3:5])
	ax_set.append(ax1)
	ax1.scatter(BE_DFT, BE_ANI, marker="o", color=color_other[0], s=markersize_value)
	ax1.scatter(BE_DFT, BE_PM7, marker="o", color=color_other[1], s=markersize_value)
	ax1.xaxis.set_major_locator(x_major_locator)
	ax1.yaxis.set_major_locator(y_major_locator)
	axis_low=min(np.min(BE_DFT), np.min(BE_ANI), np.min(BE_PM7))
	axis_high=max(np.max(BE_DFT), np.max(BE_ANI), np.max(BE_PM7))
	axis_low_about=axis_low
	axis_high_about=axis_high
	x=np.arange(axis_low_about,axis_high_about,(axis_high_about-axis_low_about)/100.0)
	y=x
	ax1.plot(x,y, color=color_other[3], linewidth=linewidth_value)
	plt.title(system_list[i-1])
	plt.xlim(axis_low_about, axis_high_about)
	plt.ylim(axis_low_about, axis_high_about) 
	ax1.text(0.015,0.93, 'RMSE_ANI-2x=%s kcal/mol'%format(RMSE_BE,'.2f'), transform=ax1.transAxes,fontsize=fontsize)
	ax1.text(0.015,0.85, 'RMSE_PM7=%s kcal/mol'%format(RMSE_BE_PM7,'.2f'), transform=ax1.transAxes,fontsize=fontsize)
	plt.xlabel('BE from DFT (kcal/mol)')
	plt.ylabel('BE from ANI-2x/PM7 (kcal/mol)')
fig1.legend(ax_set, labels=["ANI-2x", "PM7"],loc="upper right")
plt.tight_layout()	
plt.savefig('./fig2_SP_on_BH_plus_PM7.jpeg',dpi=dpi_numb)	