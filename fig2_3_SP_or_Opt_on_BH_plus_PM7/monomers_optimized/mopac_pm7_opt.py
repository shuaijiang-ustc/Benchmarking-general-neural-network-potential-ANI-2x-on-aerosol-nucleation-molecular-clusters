#import sys
import os

current_dir=os.getcwd()
a=os.listdir(current_dir)
mops=[]
for i in range(len(a)):
	if a[i][-3:]=="mop":
		mops.append(a[i])
for mop in mops:                                                                                       
	os.system('/public/scc/mopac/mopac2016/MOPAC2016.exe %s'%mop)                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                                                                                                                                                                 	                                                                                                                                                                                                                                                         	                                                                                                                                                                                                                                                           