#! /usr/bin/env python
#
# find the optimal parameters in the exchange and correlation parts of M05 functional 
#
import os,sys,string,math,time
import numpy as np
import re
import praxis
import pdb
#import scipy.optimize.optimize
#
###################################################
def prep_par(param):
#
#
# prepare the dat files with the parameters for M05
	paramx=[param[0],1.0,param[1],0.0,0.0,0.0,1.0,param[2],param[3],param[4],param[5]]
	paramcp= [1.0,param[6],param[7],0.0,0.0,param[8]]
	paramcap=[1.0,param[9],param[10],0.0,0.0]
#
# the parameters as in the original M05
#
#	 at1=	 0.08151
#	 at2=	 -0.43956
#	 at3=	 -3.22422
#	 at4=	 2.01819
#	 at5=	 8.79431
#	 at6=	 -0.00295
#	 at7=	 9.82029
#	 at8=	 -4.82351
#	 at9=	 -48.17574
#	 at10=	 3.64802
#	 at11=	 34.02248
#	 hfx =	 0.28
#	 C1 =	 3.36116E-03
#	 C2=	 4.49267E-03
#
	coefs=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
	for i in range(len(paramx)):
		coefs[i]=paramx[i]
#	coefs[12]=param[8]
# 
# coefs for revPBE
#
#	coefs[13]=0.00336116  
#	coefs[14]=0.00290129
#
	s=""
	for i in range(11):
		s+=str(coefs[i])+" "
	s+="\n"
	f=open('ParamM05X.dat','w')
	f.write(s)
#	line=f.readline()
	f.close()
#	f=open('ParamM05X.dat','r')
#	a=f.readline()
#	f.close()
#	fs=open('atta.dat','w')
#	fs.write(str(a))
#
# the parameters for the parallel part of corr in M05
#  3.77344 -26.04463 30.69913 -9.22695
#
	coefs=[0.,0.,0.,0.,0.,0.]
	for i in range(len(paramcp)):
		coefs[i]=paramcp[i]
	s=""
	for i in range(6):
		s+=str(coefs[i])+" "
	s+="\n"
	f=open('ParamM05CP.dat','w')
	f.write(s)
#	line=f.readline()
	f.close()
#	f=open('ParamM05CP.dat','r')
#	b=f.readline()
#	f.close
#	fs=open('atta.dat','w')
#	fs.write("str(b)\n")
#
# the parameters for the antiparallel part of corr in M05
# 3.78569 -14.15261 -7.46589 17.94491
#
	coefs=[0.,0.,0.,0.,0.]
	for i in range(len(paramcap)):
		coefs[i]=paramcap[i]
	s=""
	for i in range(5):
		s+=str(coefs[i])+" "
	s+="\n"
	f=open('ParamM05CAP.dat','w')
	f.write(s)
#	line=f.readline()
	f.close()
#	f=open('ParamM05CAP.dat','r')
#	c=f.readline()
#	f.close
#	fs=open('atta.dat','w')
#	fs.write(str(a))
#	fs.write(str(b))
#	fs.write(str(c)+'\n')
#	fs.close
#
	return
##################################################
def prep_run(name):
	f=open("run",'r')
	s=f.read()
	f.close()
	s=s.replace("name",name)
	f=open(prefix+"run"+name,'w')
	f.write(s)
	f.close()
	return
###################################################
def prep_inputs():
#
# prepare the inputs 
#
# herer!!!
	s0="%Mem=6GB"
	s0+="\n%NProcShared=16"
	s0+="\n#T M05/Gen test Massage SCF=(tight,maxcyc=20)"
# Guess=TCheck "
# herer!!!

#	s0+="integral=(grid=ultrafine)\n\nM05 opt"
	s0+="\n\nM05 opt"
#
	basis="\n@aug-cc-pVTZ.gbs/N\n\n"
#	count=0
#	for R in RAr:
#		count+=1
#		name="ArA_"+str(count)
#		f=open(name+".com",'w')
#		f.write("%Chk="+name+"\n")
#		f.write(s0)
#		f.write("\n\n0 1\nAr 0.0 0.0 0.0\nAr-Bq 0.0 0.0 R\n\n")
#		f.write("R="+str(R))
#		f.write("\n"+basis)
#		f.close()
#		prep_run(name)
#
# Ar2
#
#	count=0
#	for R in RAr:
#		count+=1
#		name="Ar_"+str(count)
#		f=open(name+".com",'w')
#		f.write("%Chk="+name+"\n")
#		f.write(s0)
#		f.write("\n\n0 1\nAr 0.0 0.0 0.0\nAr 0.0 0.0 R\n\n")
#		f.write("R="+str(R))
#		f.write("\n"+basis)
#		f.close()
#		prep_run(name)
#
# H2OA, H2OB, H2O 
#
	os.system('cp geoparmH2O.d geoparm.d')
	os.system('cp dimerH2O.cnf dimer.cnf')
	f=open('geoparm.d','r')
	geo=[]	
	for line in f:
		geo.append(line[0:74])	
	f.close()
#
	count=0
	for p in pH2O:
		count+=1
		ft=open('temp','w')
		ft.write(geo[p-1]+'DOR6')	
		ft.close()
		os.system('./getgeoG03.exe < temp')	
#
		name="H2OA_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoA.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="H2OB_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoB.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="H2O_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geo.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
# CH4A,CH4B,CH4 
#   
	os.system('cp geoparmCH4.d geoparm.d')
	os.system('cp dimerCH4.cnf dimer.cnf')
	f=open('geoparm.d','r')
	geo=[]	
	for line in f:
		geo.append(line[0:74])
	f.close()
#	
	count=0 
	for p in pCH4:	
		count+=1
		ft=open('temp','w')
		ft.write(geo[p-1]+'DOR6')
		ft.close()	
		os.system('./getgeoG03.exe < temp')
#			
		name="CH4A_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoA.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="CH4B_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoB.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="CH4_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geo.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
# C2H4A,C2H4B,C2H4 
#	
	os.system('cp geoparmC2H4.d geoparm.d')
	os.system('cp dimerC2H4.cnf dimer.cnf')
	f=open('geoparm.d','r')
	geo=[]
	for line in f:
		geo.append(line[0:74])
	f.close()
#	
	count=0
	for p in pC2H4:
		count+=1
		ft=open('temp','w')
		ft.write(geo[p-1]+'DOR6')
		ft.close()
		os.system('./getgeoG03.exe < temp')
#			
		name="C2H4A_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write("Density=Checkpoint\n")
#		f.write("%OldChk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoA.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="C2H4B_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write("Density=Checkpoint\n")
#		f.write("%OldChk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoB.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="C2H4_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write("Density=Checkpoint\n")
#		f.write("%OldChk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geo.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
# C2H2ClFA, C2H2ClFB, C2H2ClF
#
	os.system('cp geoparmC2H2ClF.d geoparm.d')
	os.system('cp dimerC2H2ClF.cnf dimer.cnf')
	f=open('geoparm.d','r')
	geo=[]
	for line in f:
		geo.append(line[0:74])
	f.close()
#
	count=0
	for p in pC2H2ClF:
		count+=1
		ft=open('temp','w')
		ft.write(geo[p-1]+'DOR6')
		ft.close()
		os.system('./getgeoG03.exe < temp')
#
		name="C2H2ClFA_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoA.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="C2H2ClFB_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoB.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="C2H2ClF_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geo.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
# NH3_ClFA, NH3_ClFB, NH3_ClF
#
	os.system('cp geoparmNH3_ClF.d geoparm.d')
	os.system('cp dimerNH3_ClF.cnf dimer.cnf')
	f=open('geoparm.d','r')
	geo=[]
	for line in f:
		geo.append(line[0:74])
	f.close()
#
	count=0
	for p in pNH3_ClF:
		count+=1
		ft=open('temp','w')
		ft.write(geo[p-1]+'DOR6')
		ft.close()
		os.system('./getgeoG03.exe < temp')
#
		name="NH3_ClFA_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoA.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="NH3_ClFB_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoB.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="NH3_ClF_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geo.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
# C2H4_ArA, C2H4_ArB, C2H4_Ar
#
	os.system('cp geoparmC2H4_Ar.d geoparm.d')
	os.system('cp dimerC2H4_Ar.cnf dimer.cnf')
	f=open('geoparm.d','r')
	geo=[]
	for line in f:
		geo.append(line[0:74])
	f.close()
#
	count=0
	for p in pC2H4_Ar:
		count+=1
		ft=open('temp','w')
		ft.write(geo[p-1]+'DOR6')
		ft.close()
		os.system('./getgeoG03.exe < temp')
#
		name="C2H4_ArA_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoA.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="C2H4_ArB_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoB.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="C2H4_Ar_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geo.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
# C2H4_F2A, C2H4_F2B, C2H4_F2 
#	
#	os.system('cp geoparmC2H4_F2.d geoparm.d')
#	os.system('cp dimerC2H4_F2.cnf dimer.cnf')
#	f=open('geoparm.d','r')
#	geo=[]
#	for line in f:
#		geo.append(line[0:74])
#	f.close()
#	
#	count=0
#	for p in pC2H4_F2:
#		count+=1
#		ft=open('temp','w')
#		ft.write(geo[p-1]+'DOR6')
#		ft.close()
#		os.system('./getgeoG03.exe < temp')
#			
#		name="C2H4_F2A_"+str(count)
#		f=open(name+".com",'w')
#		f.write("%Chk="+name+"\n")
#		f.write(s0)
#		f.write("\n\n0 1\n")
#		fgeo=open('geoA.d','r')
#		s=fgeo.read()
#		fgeo.close()
#		f.write(s)
#		f.write(basis)
#		f.close()
#		prep_run(name)
#
#		name="C2H4_F2B_"+str(count)
#		f=open(name+".com",'w')
#		f.write("%Chk="+name+"\n")
#		f.write(s0)
#		f.write("\n\n0 1\n")
#		fgeo=open('geoB.d','r')
#		s=fgeo.read()
#		fgeo.close()
#		f.write(s)
#		f.write(basis)
#		f.close()
#		prep_run(name)
#
#		name="C2H4_F2_"+str(count)
#		f=open(name+".com",'w')
#		f.write("%Chk="+name+"\n")
#		f.write(s0)
#		f.write("\n\n0 1\n")
#		fgeo=open('geo.d','r')
#		s=fgeo.read()
#		fgeo.close()
#		f.write(s)
#		f.write(basis)
#		f.close()
#		prep_run(name)
#
# NH3A, NH3B, NH3
#
	os.system('cp geoparmNH3.d geoparm.d')
	os.system('cp dimerNH3.cnf dimer.cnf')
	f=open('geoparm.d','r')
	geo=[]
	for line in f:
		geo.append(line[0:74])
	f.close()
#
	count=0
	for p in pNH3:
		count+=1
		ft=open('temp','w')
		ft.write(geo[p-1]+'DOR6')
		ft.close()
		os.system('./getgeoG03.exe < temp')
#
		name="NH3A_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoA.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="NH3B_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoB.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="NH3_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geo.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
# ArHFA, ArHFB, ArHF
#
	os.system('cp geoparmArHF.d geoparm.d')
	os.system('cp dimerArHF.cnf dimer.cnf')
	f=open('geoparm.d','r')
	geo=[]
	for line in f:
		geo.append(line[0:74])
	f.close()
#
	count=0
	for p in pArHF:
		count+=1
		ft=open('temp','w')
		ft.write(geo[p-1]+'DOR6')
		ft.close()
		os.system('./getgeoG03.exe < temp')
#
		name="ArHFA_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoA.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="ArHFB_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoB.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="ArHF_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geo.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
# N2_HFA, N2_HFB, N2_HF
#
	os.system('cp geoparmN2_HF.d geoparm.d')
	os.system('cp dimerN2_HF.cnf dimer.cnf')
	f=open('geoparm.d','r')
	geo=[]
	for line in f:
		geo.append(line[0:74])
	f.close()
#
	count=0
	for p in pN2_HF:
		count+=1
		ft=open('temp','w')
		ft.write(geo[p-1]+'DOR6')
		ft.close()
		os.system('./getgeoG03.exe < temp')
#
		name="N2_HFA_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoA.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="N2_HFB_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoB.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="N2_HF_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geo.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
		prep_run(name)
#
# change the basis set to aDZ
#
	basis="\n@aug-cc-pVDZ.gbs/N\n\n"
#
# H2ODZA, H2ODZB, H2ODZ 
#
	os.system('cp geoparmH2OZ.d geoparm.d')
	os.system('cp dimerH2OZ.cnf dimer.cnf')
	f=open('geoparm.d','r')
	geo=[]
	for line in f:
		geo.append(line[0:74])
	f.close()
#
	count=0
	for p in pH2ODZ:
		count+=1
		ft=open('temp','w')
		ft.write(geo[p-1]+'DOR6')
		ft.close()
		os.system('./getgeoG03.exe < temp')
#
		name="H2ODZA_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoA.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="H2ODZB_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geoB.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
		name="H2ODZ_"+str(count)
		f=open(name+".com",'w')
		f.write("%Chk="+name+"\n")
		f.write(s0)
		f.write("\n\n0 1\n")
		fgeo=open('geo.d','r')
		s=fgeo.read()
		fgeo.close()
		f.write(s)
		f.write(basis)
		f.close()
		prep_run(name)
#
	os.system("rm temp")
#
	return
###################################################
def en_int(monomerA,monomerB,dimer):
#
	conv=627.509474
# get the energy of a monomer in DCBS 
#
	Conv1=0
	f=open(monomerA,'r')
	for line in f:
		if line[:10]==" SCF Done:":
			Conv1=1
			for l in range(len(line)) :
				if line[l]=="=":
					ini=l+1
				if line[l]=="A":
					fin=l-2
			en_monoA=float(line[ini:fin])
#			print "en_monoA", en_monoA
	f.close()
	if Conv1==0:
		en_monoA=0.0
		print("No convergence in"),monomerA
#		print "No convergence in",monomerA
#
	Conv1=0
	f=open(monomerB,'r')
	for line in f:
		if line[:10]==" SCF Done:":
			Conv1=1
			for l in range(len(line)) :
				if line[l]=="=":
					ini=l+1
				if line[l]=="A":
					fin=l-2
			en_monoB=float(line[ini:fin])
#			print "en_monoB", en_monoB
	f.close()
	if Conv1==0:
		en_monoB=0.0
		print("No convergence in"),monomerB
#		print "No convergence in",monomerB
#
#
# obtain the energy of dimer
#
	Conv2=0
	f=open(dimer,'r')
	for line in f:
		if line[:10]==" SCF Done:":
			Conv2=1
			for l in range(len(line)) :
				if line[l]=="=":
					ini=l+1
				if line[l]=="A":
					fin=l-2
			en_dim=float(line[ini:fin])
	f.close()
	if Conv2==0:
		en_dim=0.0
		print("No convergence in"),dimer
#		print "No convergence in",dimer
#
#	Dispersion part:
#
#	file='/data/szalewicgrp/atta/dldf-praxis-opt/Disp/'+monoA
#	f=open('/data/szalewicgrp/atta/dldf-test-praxis-opt/monomers/'+monoA,'r')
#	lineA=f.read()
#	f3=open('fileA.dat','w')
#	f3.write(str(lineA))
#	f3.close()
#	f.close()
#	file='/data/szalewicgrp/atta/dldf-praxis-opt/Disp/'+monoB
#	f=open('/data/szalewicgrp/atta/dldf-test-praxis-opt/monomers/'+monoB,'r')
#	lineB=f.read()
#	f4=open('fileB.dat','w')
#	f4.write(lineB)
#	f4.close()
#	f.close()
#	os.system('./disp.exe')
#	with open('aaa','r') as f:
#		d_e=f.readline()
#	disp_en=float(d_e)
#	fd=open('/data/szalewicgrp/atta/dldf-test-praxis-opt/monomers/disp.data','a')
#	fd.write(monoA+'\n')
#	fd.write(monoB+'\n')
#	fd.write(str(disp_en)+'\n')
#	fd.close()
	dlDF=round((en_dim-en_monoA-en_monoB)*conv,8)
	return dlDF
###################################################
def exist(name):
	f=os.popen("ls "+name+" 2>sss")
	s=f.read()
	f.close()
	os.system("rm sss")
	return s==name+"\n"
###################################################
def get_queue():
	listq=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
	qno=-1
	while qno==-1:
		for q in listq:
			f=os.popen("qstat -f")
			s=f.read()
			f.close()
			ind=s.find("all.q@compute-0-"+str(q))
			used=float(s[ind+37:ind+38])
			f=os.popen("qstat -q all.q@compute-0-"+str(q)+".local | grep '	qw  '")
			waiting=f.read()
			f.close()
			f=os.popen("qstat -q all.q@compute-0-"+str(q)+".local | grep '	t  '")
			t=f.read()
			f.close()
#
			if q==14:
				used+=1
#
			if used<2 and len(waiting)==0 and len(t)==0:
				qno=q
				break
	return qno
###################################################
def jobs_running(m):
	run=0
	while run==0:
		f=os.popen("squeue -u connornm | grep 'connornm'")
		count=0
		for line in f:
			count+=1
		f.close()
		if count<=m:
			run=1
	return 
###################################################
def get_err(param,n):
	global x
	x= param
#	param=x[0]
#	conv=627.509474
	n=11
#get_err(param,names,ERef,RRef,no_pts)
#	conv=3.157733E5
#
	prep_par(param)
#
# prepare a list of jobs
# here range(3) means 3 molecules, change for more molecules
	list=[]
	for m in range(10):
		count=0
		for i in range(no_pts[m]):
			count+=1
			name="run"+names[m]+"A_"+str(count)
			list.append(prefix+name)
			name="run"+names[m]+"B_"+str(count)
			if exist(prefix+name):
				list.append(prefix+name)
			name="run"+names[m]+"_"+str(count)
			list.append(prefix+name)
#
	for job in list:
#		jobs_running(5)
#		print job
#		q=get_queue()
#		os.system("qsub -q all.q@compute-0-"+str(q)+".local "+job+" >sss")
		os.system("sbatch "+job+" >sss")
	os.system("rm sss")
#
	running=1
	while running:
		f=os.popen("squeue -u connornm")
		s=f.read()
		f.close()
		ind=s.find(prefix+"run")
		if ind==-1:
			running=0
# collect the results
#
	eDFT=[]
	molec=[]
#
# add the errors from the interaction energies
#
	for n in range(len(names)):
		count=0	
		for i in range(no_pts[n]):
			count+=1
			monomerA=names[n]+"A_"+str(count)+".log"
			monomerB=names[n]+"B_"+str(count)+".log"
#			monoA=names[n]+"_"+str(count)+"_monomerA.com"
#			monoB=names[n]+"_"+str(count)+"_monomerB.com"
			if not exist(monomerB):
				monomerB=monomerA
			dimer= names[n]+"_"+str(count)+".log"
			eDFT.append(en_int(monomerA,monomerB,dimer))
			molec.append(names[n])
#
	print("Dispersionless energies")
#	print "Dispersionless energies"
	print("{0:7}     {1:14}  {2:14}   {3:14}   {4:14}".format('Molecule','R[a.u.]','E_DFT[K]','E_Ref[K]','Err'))
#	print "R[a.u.], E_DFT[K], E_Ref[K], Err"
	err=0
	error=0
	fe=open('energy.dat','a')
	fe.write("{0:10}  {1:14}  {2:14}   {3:14}   {4:14}   {5:14}\n".format('Molecule','R[a.u.]','E_DFT[Kcal/mol]','E_Ref[Kcal/mol]','Abs_Err[Kcal/mol]','SQ_Rel_Err'))
	for i in range(len(ERef)):
#		print(eDFT[i])
		e_ref=ERef[i]
		e_dft=round(eDFT[i],8)
		R=RRef[i]
		if (e_ref!=0):
#
#			err+=((e_dft-e_ref)/min(abs(e_ref),abs(e_dft)))**2
			err+=((e_dft-e_ref)/e_ref)**2
			error+=((e_dft-e_ref)/e_ref)**2
		if (e_ref!=0):
#			print( molec[i],R,e_dft,e_ref,((e_dft-e_ref)/min(abs(e_ref),abs(e_dft)))**2)
#			print( molec[i],R,e_dft,e_ref,((e_dft-e_ref)/e_ref)**2)
			print("{0:7} {1:14} {2:14} {3:14} {4:14}".format(molec[i],round(R,5),e_dft,e_ref,round(((e_dft-e_ref)/e_ref)**2,8)))
			fe.write("{0:7} {1:14} {2:14} {3:18} {4:18} {5:18}\n".format(molec[i],round(R,5),e_dft,e_ref,round(abs(e_dft-e_ref),6),round(((e_dft-e_ref)/e_ref)**2,8)))
		else:
			print( molec[i],R,e_dft,e_ref)
#	fe.write("*****************************NEW Iteration*****************************" + '\n')
#	fe.close()
#
# add the errors from the dipole moments
#
	print("Dipole Moments")
	print("Dipole_DFT, Dipole_Ref, Err")
	fe.write("{0:10}   {1:14} {2:14} {3:14}\n".format('Molecule','Dipole_DFT', 'Dipole_Ref', 'Err'))
#
	weightDip=100.
#
# from H2OA_3	

# S in the text with the whole file in it
#
	Conv1=0
	f=open("H2OA_3.log","r")
	s=f.read()
	f.close()
	f=open("H2OA_3.log","r")
	for line in f:
		if line[:20]==" Convergence failure":
			Conv1=1
#		if line[:10]==" SCF Done:":
#			Conv1=2
	f.close()
	if Conv1==1:
		err+=10
#       Don't understand this part of the code -- Connor vvvvvvv
#	if Conv1==2:
	else:
	#  Find X from log file from Gaussian
		ind=s.find("    X=")
		dip_dft=float(s[ind+19:ind+26])
		# Find difference in measured dipole and dipole reference value 
		err+=((dip_dft-DipRef[0])/DipRef[0])**2*weightDip
		error+=((dip_dft-DipRef[0])/DipRef[0])**2
		print( "H2O_X",dip_dft,DipRef[0],((dip_dft-DipRef[0])/DipRef[0])**2*weightDip)
		fe.write("{0:7} {1:14} {2:14} {3:14}\n".format('H2O_X',dip_dft,DipRef[0],round(((dip_dft-DipRef[0])/DipRef[0])**2*weightDip,8)))
#	
		ind=s.find("Z=")
		dip_dft=float(s[ind+13:ind+22])
		err+=((dip_dft-DipRef[1])/DipRef[1])**2*weightDip
		error+=((dip_dft-DipRef[1])/DipRef[1])**2
		print( "H2O_Z",dip_dft,DipRef[1],((dip_dft-DipRef[1])/DipRef[1])**2*weightDip)
		fe.write("{0:7} {1:14} {2:14} {3:14}\n".format('H2O_Z',dip_dft,DipRef[1],round(((dip_dft-DipRef[1])/DipRef[1])**2*weightDip,8)))
#
#	else:
#		print ("No convergence in H2OA_3.log")
#		err+=100000
#
	Conv1=0
	f=open("NH3A_3.log","r")
	s=f.read()
	f.close()
#
	f=open("NH3A_3.log","r")
	for line in f:
		if line[:20]==" Convergence failure":
			Conv1=1
#		if line[:10]==" SCF Done:":
#			Conv1=1
#
	f.close()
	if Conv1==1:
		err+=10
	else:
		ind=s.find("Y=")
		dip_dft=float(s[ind+13:ind+24])
		err+=((dip_dft-DipRef[2])/DipRef[2])**2*weightDip
		error+=((dip_dft-DipRef[2])/DipRef[2])**2
		print( "NH3_Y",dip_dft,DipRef[2],((dip_dft-DipRef[2])/DipRef[2])**2*weightDip)
		fe.write("{0:7} {1:14} {2:14} {3:14}\n".format('NH3_Y',dip_dft,DipRef[2],round(((dip_dft-DipRef[2])/DipRef[2])**2*weightDip,8)))
#	
		ind=s.find("Z=")
		dip_dft=float(s[ind+13:ind+22])
		err+=((dip_dft-DipRef[3])/DipRef[3])**2*weightDip
		error+=((dip_dft-DipRef[3])/DipRef[3])**2
		print( "NH3_Z",dip_dft,DipRef[3],((dip_dft-DipRef[3])/DipRef[3])**2*weightDip)
		fe.write("{0:7} {1:14} {2:14} {3:14}\n".format('NH3_Z',dip_dft,DipRef[3],round(((dip_dft-DipRef[3])/DipRef[3])**2*weightDip,8)))
#	else:
#		print ("No convergence in H2OA_3.log")
#		err+=100000
#	fe.write("*****************************New Loop*****************************" + '\n')
#	fe.close()
	error=math.sqrt(error/33)*(100)
	err=math.sqrt(err/33)*(100)
	fe.write("{0:10} {1:10} {2:10} {3:10}\n".format("Error_Weight%=",err,"Error_NO_Weight%=",error))
	fe.write("*****************************New Loop*****************************" + '\n')
	fe.close()
#	get_err=err
#	time.sleep(30)
	print(error)
	print(param,'Percent Error:',err)
	f=open('param.dat','a')
#	f.write(str(error) + '\n')
	f.write("{0:14} {1:14} {2:14} {3:14} {4:14} {5:14} {6:14}\n".format(str(param[0]),str(param[1]),str(param[2]),str(param[3]),str(param[4]),\
	str(param[5]),str(param[6]),str(param[7]),str(param[8]),str(param[9]),str(param[10]),'% Error_weight:',round(err,8),'% Error:',round(error,8)))
	f.close()
#	for i in param:
#		print(i)
#	
	if err < 11:
		sys.exit(0)
	return err
####################################################
def prep_data():
#
# conversion factor from Hartree to K
#       conv_fac=3.157733E5
	conv_fac=1
	print( "R[a.u.]  Eint[Kcal/mol]")
#
# All R in [A], energies in [Hartree] from SAPT(DFT)_dispersionless
#
# H2O
	R_All=[2.5, 3.1, 4.0]
	ERef_All=[5.297988310, -3.023071530, -1.766518960]
	if len(pH2O)>0:
		print("H2O")
	for i in pH2O:
		RH2O.append(R_All[i-1])
		RRef.append(R_All[i-1])
		ERefH2O.append(ERef_All[i-1])
		ERef.append(ERef_All[i-1])
		print(R_All[i-1]/0.529177208,ERef_All[i-1]*conv_fac)
#
# for aDZ basis set
#
# H2ODZ
#
#	R_All=[2.5, 3.1, 4.0]
#	ERef_All=[5.457619510, -3.035855110, -1.776455210]
	R_All=[3.1, 4.0]
	ERef_All=[-3.035855110, -1.776455210]
	if len(pH2ODZ)>0:
		print("H2ODZ")
	for i in pH2ODZ:
		RH2ODZ.append(R_All[i-1])
		RRef.append(R_All[i-1])
		ERefH2ODZ.append(ERef_All[i-1])
		ERef.append(ERef_All[i-1])
		print(R_All[i-1]/0.529177208,ERef_All[i-1]*conv_fac)
#
#
# Reference dipole moments in [Debye] from CCSD(T) aDZ 
#
# for H2OA_3
	DipRef.append(-1.61256588)
	DipRef.append(0.89657227)
#
# for NH3A_3
	DipRef.append(1.41009079)
	DipRef.append(-0.57028608)
#
	return
#########################################
#     BEGINNING OF THE MAIN PART	#
#########################################
#
# input the initial parameters 
#
# set prefix for job names
prefix="S"
#
# the parameters come from the Genetic Algorithm Opt 
# he 9th parameter corresponds to xhf!
#restart("mn_opt_chimera.py")
#param=[-7.83343246615, 15.26130438, 1.38517538]
#param=[-7.84391049533, 15.26130438, 1.36517538]
#param=[-7.83343247, 15.26130438, 3.38517538]
#param=[-7.31481632, 15.57423523, 5.35533886]    #76%
#param=[-1.16011578, 19.0666469, 2.20063832]
#param=[1.851851, 17.33, 0.163]
param=[1.18328602E-03, 4.8827323, -0.163757125, -0.188002773, -0.449060892, -0.00823585719, 5.9515308, \
	-11.1602877, 0.614413, -2.59608972, 2.22337926]
#param=[1.38303837e-03, 4.09322728e+00, -3.32922779e-01, -2.44905609e-01, -3.91247752e-01, -4.95542667e-02, \
#	4.56238314e+00, -1.08236169e+01, 5.72471781e-01, -2.96133418e+00, 2.04022220e+00]
#param=[1.18328602E-03, 3.8827323E+00, -1.63757125E-01, -1.88002773E-01, -3.49060892E-01, -7.23585719E-03, \
#	3.95153080E+00, -1.01602877E+01, 0.56, -2.59608972E+00, 2.22337926E+00]
#param=[ 0.00891661528259, 3.93680079715, -0.10684221087, -0.35860210887, -0.524278554611, -0.56236347868, \
#	5.52213013587, -9.58968836413, 0.775493832148, -3.06668905587, 1.75277992413]
#param=[ 0.05281973, 0.54037798, -0.69588941, -0.19629029, -2.62045642, -4.95992745, \
#	-1.91062223, 1.21189836, 0.12786006, 4.23421, -8.24132]
#param=[ 1.81917515, 0.0864756, -3.6452664000000001, -1.65665967,  \
#	0.095848550000000005, -2.66675168, -0.35652831000000001, -2.82758398, 5.5, -10.2]
#
# reference points pxxx=[] means exluding xxx form the optimization set 
#
pC2H4=[1,2,3]
pC2H2ClF=[1,2,3]
pNH3_ClF=[1,2,3]
pC2H4_Ar=[1,2]
pH2O=[1,2,3]
pCH4=[1,2,3]
pNH3=[1,2,3]
pArHF=[1,2,3]
pN2_HF=[1,2,3]
# for aDZ basis set
pH2ODZ=[1,2]
#pHe=[]
#pHe=[2,6,10]
#pAr=[8,12,17]
#pAr=[]
#pH2O=[]
#pH2O=[2,4,6]
#pH2O=[]
#pCH4=[]
#pCH4=[1,4,6]
#pCH4=[]
#pC2H4=[]
#pC2H4=[2,4,7]
#pC2H4=[]
#pC2H4_F2=[3,5,7]
#pC2H4_F2=[]
#pNH3=[]
#pNH3=[2,3,6]
#pArHF=[]
#pArHF=[2,3,5]
#pN2_HF=[]
#pN2_HF=[2,4,7]
# for aDZ basis set
#pH2ODZ=[2,4,6]
names=["C2H4","C2H2ClF","NH3_ClF","C2H4_Ar","H2O","CH4","NH3","ArHF","N2_HF","H2ODZ"]
no_pts=[len(pC2H4),len(pC2H2ClF),len(pNH3_ClF),len(pC2H4_Ar),len(pH2O),len(pCH4),len(pNH3), \
	len(pArHF),len(pN2_HF),len(pH2ODZ)]
#########################################
#
n=11
ERef=[]
RRef=[]
RC2H4=[]
ERefC2H4=[]
RC2H2ClF=[]
ERefC2H2ClF=[]
RNH3_ClF=[]
ERefNH3_ClF=[]
RC2H4_Ar=[]
ERefC2H4_Ar=[]
RH2O=[]
ERefH2O=[]
RCH4=[]
ERefCH4=[]
RNH3=[]
ERefNH3=[]
RArHF=[]
ERefArHF=[]
RN2_HF=[]
ERefN2_HF=[]
RH2ODZ=[]
ERefH2ODZ=[]
#
DipRef=[]
#
prep_data()
prep_inputs()
#
os.system("rm *.chk")
#
# optimize the parametres 
#
x=praxis.praxis(0.00001,0.00001,11,4,param,get_err)
#time.sleep(30)
print(x)
param=x
#f=open('res1.data','w')
#f.write(str(x[0]))
#f.close()
#sys.exit()
#
# print the interaction energies for all distances
#
os.system("rm *.chk")
#pHe=range(1,15)
#pAr=range(1,24)
pC2H4=range(1,4)
pC2H2ClF=range(1,4)
pNH3_ClF=range(1,4)
pC2H4_Ar=range(1,3)
pH2O=range(1,4)
pCH4=range(1,4)
pNH3=range(1,4)
pArHF=range(1,4)
pN2_HF=range(1,4)
# for aDZ basis set
pH2ODZ=range(1,3)
#
#pC2H4=[]
#pC2H2ClF=[]
#pNH3_ClF=[]
#pC2H4_Ar=[]
#pH2O=[]
#pCH4=[]
#pNH3=[]
#pArHF=[]
#pN2_HF=[]
# for aDZ basis set
#pH2ODZ=[]
#
names=["C2H4","C2H2ClF","NH3_ClF","C2H4_Ar","H2O","CH4","NH3","ArHF","N2_HF","H2ODZ"]
no_pts=[len(pC2H4),len(pC2H2ClF),len(pNH3_ClF),len(pC2H4_Ar),len(pH2O),len(pCH4),len(pNH3), \
	len(pArHF),len(pN2_HF),len(pH2ODZ)]
#
ERef=[]
RRef=[]
RC2H4=[]
ERefC2H4=[]
RC2H2ClF=[]
ERefC2H2ClF=[]
RNH3_ClF=[]
ERefNH3_ClF=[]
RC2H4_Ar=[]
ERefC2H4_Ar=[]
RH2O=[]
ERefH2O=[]
RCH4=[]
ERefCH4=[]
RNH3=[]
ERefNH3=[]
RArHF=[]
ERefArHF=[]
RN2_HF=[]
ERefN2_HF=[]
RH2ODZ=[]
ERefH2ODZ=[]
#
DipRef=[]
#
prep_data()
prep_inputs()
#
print( "FINAL RESULTS")
err=get_err(param,n)
#
#sys.exit()
