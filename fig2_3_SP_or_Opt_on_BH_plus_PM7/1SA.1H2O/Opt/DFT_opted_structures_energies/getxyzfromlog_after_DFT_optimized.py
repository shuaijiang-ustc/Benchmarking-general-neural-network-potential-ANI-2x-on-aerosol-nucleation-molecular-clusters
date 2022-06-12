import os
import sys

current_dir=os.getcwd()
a=os.listdir(current_dir)
log=[]
for i in range(len(a)):
	if a[i][-3:]=="log":
		log.append(a[i])

ab=[]
abc=[]
for line in log:
	outname=line
	nor=os.popen('cat %s | grep "Normal termination"'%outname).readlines()
	if len(nor)==2:
		a=os.popen('cat %s | grep "NAtoms="'%outname).readlines()
		number=int(a[0].split()[1])
		b=os.popen('cat -n %s | grep "Standard orientation:"'%outname).readlines()
		bb=os.popen('cat -n %s | grep "Symbolic Z-matrix:"'%outname).readlines()
		bbb=os.popen('cat -n %s | grep "SCF Done:"'%outname).readlines()
		firstnumber=int(b[-1].split()[0])
		elementsymbonumber=int(bb[0].split()[0])
		energy=bbb[-1].split()[5]
		ab.append(float(energy))
		abc.append("Geom_"+line[:-4])
		n1=firstnumber+5
		n2=firstnumber+number+4
		n11=elementsymbonumber+2
		n22=elementsymbonumber+number+1
		d=os.popen('cat %s | sed -n "%s,%sp"'%(outname,n1,n2)).readlines()
		dd=os.popen('cat %s | sed -n "%s,%sp"'%(outname,n11,n22)).readlines()
		geometry=""
		geometry+=str(number)+"\n"
		geometry+=str(energy)+"\n"
		for i in range(len(d)):
			geometry+=dd[i].split()[0]+"  "+d[i].split()[3]+"  "+d[i].split()[4]+"  "+d[i].split()[5]+"\n"
			xyzname='DFT_optimized_'+line[:-4]+".xyz"
			xyz=open(xyzname,"w")
			xyz.write(geometry)
			xyz.close()
	else:
		pass

without_sort=zip(ab,abc)
without_sort.sort()

relatenergy=[]
for i in range(len(without_sort)):
	a=(without_sort[i][0]-without_sort[0][0])*627.5094
	relatenergy.append(a)

energystring=""
ss=0
for i in range(len(without_sort)):
	energystring+=str(ss)+"  "+str(without_sort[i][1])+"  "+str(without_sort[i][0])+"  "+str(relatenergy[i])+"\n"
	ss=ss+1
	
energyfile=open("energy_sorted.txt","w")
energyfile.write('order'+'\t'+'isomer_name'+'\t'+'energy_in_Hatree'+'\t'+'energy_difference_in_kcal/mol'+'\n')
energyfile.write(energystring)
energyfile.close()
