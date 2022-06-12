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
x_major_locator=MultipleLocator(15)
y_major_locator=MultipleLocator(15)
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
dir_list=['SA_W\\bowl_0.05', 'SA_Am\\bowl_0.05', 'SA_DMA\\bowl_0.05', 'SA_SA\\bowl_0.05', 'SA_SUA\\bowl_0.05']
system_list=['(SA)$_{1}$(W)$_{1}$', '(SA)$_{1}$(Am)$_{1}$','(SA)$_{1}$(DMA)$_{1}$','(SA)$_{2}$','(SA)$_{1}$(SUA)$_{1}$' ]

i=0
ax_set=[]
fig1=plt.figure(num=1, figsize=(8.5, 5.5), dpi=600, facecolor='w', edgecolor='k')
gs = gridspec.GridSpec(2, 6)
for dirt in dir_list:
	i+=1
	target_path=os.path.join(cur_dir, dirt)
	target_path=os.path.join(target_path, 'wB97X_6-31G_star')
	target_file=os.path.join(target_path, 'ANI_DFT_RE.txt')
	f=open(target_file,'r')
	f.readline()
	RE_DFT=[]
	RE_ANI=[]
	for line in f.readlines():
		RE_DFT.append(float(line.split()[0]))
		RE_ANI.append(float(line.split()[1]))
	f.close()
	RMSE_RE = metrics.mean_squared_error(RE_ANI, RE_DFT)**0.5
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
	ax1.scatter(RE_DFT, RE_ANI, marker="o", color=color_other[0], s=markersize_value)
	ax1.xaxis.set_major_locator(x_major_locator)
	ax1.yaxis.set_major_locator(y_major_locator)
	axis_low=min(np.min(RE_DFT), np.min(RE_ANI))
	axis_high=max(np.max(RE_DFT), np.max(RE_ANI))
	axis_low_about=axis_low
	axis_high_about=axis_high
	x=np.arange(axis_low_about,axis_high_about,(axis_high_about-axis_low_about)/100.0)
	y=x
	ax1.plot(x,y, color=color_other[3], linewidth=linewidth_value)
	plt.title(system_list[i-1])
	plt.xlim(axis_low_about, axis_high_about)
	plt.ylim(axis_low_about, axis_high_about) 
	ax1.text(0.015,0.93, 'RMSE_ANI-2x=%s kcal/mol'%format(RMSE_RE,'.2f'), transform=ax1.transAxes,fontsize=fontsize)
	plt.xlabel('RE from DFT (kcal/mol)')
	plt.ylabel('RE from ANI-2x (kcal/mol)')
	
fig1.legend(ax_set, labels=["ANI-2x"],loc="upper right")
plt.tight_layout()	
plt.savefig('./fig6_SP_on_MetaMD.jpeg',dpi=dpi_numb)	