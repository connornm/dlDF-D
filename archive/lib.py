from subprocess import check_output
import os, pdb, re


#	wait - Wait until all jobs are done	
#########################################################
def wait(systems):
	for system in systems:
		while bool(int(check_output(['bin/number_of_runs', system]))):
			pass

#	run - Run systems using specific software
#########################################################
def run(systems, software):
	for system in systems:
		os.system('./run'+software+' ' + system)

#	energy_ - Get array of energies for list of 
#	systems for a specific software keyword	
#########################################################
def energy(systems, software):
	E = {}
	for system in systems:
		try:
			E[system]=float(check_output(['bin/energy_from_'+software, system]))
		except:
			print('Could not get energy from '+system)
	return E

#	write - write a string to the output file
#########################################################
def write(instr):
	os.system("echo '"+instr+"' >> output")


#	mon_dim_energy - Get dictionary of monomer and 
#	dimer energies from Gauss	
#########################################################
#	Input:
#	systems - Name of the systems
#########################################################
#	Output:
#	E_mult - Dictionary of energies 
#########################################################
def mon_dim_energy(systems):
	E_mult = {}
	for system in systems:
		E_mult[system+'_D'] = get_energy(system+'_D')
		E_mult[system+'_A'] = get_energy(system+'_A')
		E_mult[system+'_B'] = get_energy(system+'_B')
	return E_mult

#	mon_dim_run - Run set of monomer and dimer jobs
#########################################################
#	Input:
#	system - The name of the system
#########################################################
def mon_dim_run(systems):
	for system in systems:
		os.system('./run ' + system + '_D')
		os.system('./run ' + system + '_A')
		os.system('./run ' + system + '_B')

#	save - Saves an array to an existing file
#########################################################
def save(array, filename):
	fstr = ''
	for i in range(len(array)):
		fstr += str(array[i]) + ' '
	fstr += '\n'
	f = open(filename, 'w')
	f.seek(0)
	f.write(fstr)
	f.close()	
	
#	load - Loads an array from a file
#########################################################
def load(filename):
	array = []
	f = open(filename, 'r')
	f.seek(0)
	fstr=f.read()
	f.close()
	fvect = fstr.split(' ') 
	for i in range(len(fvect)):
		try:
			array += [float(fvect[i])]
		except:
			pass	
	return array

#	read - Read arrays for systems
#########################################################
def read(systems):
	data = {}	
	for system in systems:
		temp = load('data/ener/'+system+'.ener')
		system_data = {
			# Intermolecular distance
			#'R': bench_data[0], 
			# Total energy
			'E': temp[-10],
			# Monomer A energy
			#'E_A': bench_data[-2],
			# Monomer B energy
			#'E_B': bench_data[-1]
		}
		data[system] = system_data
	return data

#	prep_dispersion
#########################################################
#	Input:
#	systems - systems to create .cnf files for from 
#	.com files
#########################################################
def prep_dispersion(systems):
	for system in systems:
		os.system('bin/com2cnf ' + system)
	
#	disp_energy
#########################################################
#	Input:
#	system - system to obtain dispersion energy from
#	gets from a .cnf file
#########################################################
#	Output:	
#	denergy - dispersion energy of system
#########################################################
def disp_energy(system):
	return float(check_output(['./disp4dldf.py','data/cnf/'+system+'.cnf', '2']))	

#	mult_disp
#########################################################
#	Input:
#	systems - list of systems to get dispersion 
#	energy from
#########################################################
#	Output:
#	eDISP - dictionary of dispersion energies with 
#	systems as the keys
#########################################################
def mult_disp(systems): 
	eDISP = {}
	for system in systems:
		eDISP[system] = disp_energy(system)
	return eDISP

#	prep_dim_mon
#########################################################
#	Input:
#	systems - list of systems to prepare dimer and 
#	monomer .com files
#########################################################
def prep_dim_mon(systems):
	for system in systems:
		os.system('bin/prep_dim_mon '+system)

#	prep_gauss
#########################################################
#	Input:
#	systems - list of systems to prepare .com files
#	from .xyz files
#########################################################
def prep_gauss(systems):
	for system in systems:
		os.system('bin/xyz2com '+system)	
