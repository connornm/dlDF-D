from subprocess import check_output
import os, pdb, re


#	get_energy - Reads energy from Gauss file
#########################################################
#	Input:
#	system - The name of the system 
#########################################################
#	Output:	
#	energy - The energy from Gaussian 
#########################################################
def get_energy(system):
	try:
		return float(check_output(['bin/energy_from_gaussian', system]))
	except:
		print('Could not get energy from system: '+system)
		return 0

#	run - Submit Gauss job for system 
#########################################################
#	Input:
#	system - The name of the system
#########################################################
def run(system):
	os.system('./run ' + system)

#	wait - Wait until a Gaussian job is done	
#########################################################
#	Input:
#	system - The name of the system 
#########################################################
def wait(system):
	while bool(int(check_output(['bin/number_of_runs', system]))):
		pass

#	read_param - Reads parameters from param file
#########################################################
#	Input:
#	func - The name of the functional to paramereterize
#	init - Boolean where init=True indicates reading
#	parameters from the initial directory 
#########################################################
#	Output: 
#	params - Parameters of the M05X	
#########################################################
def read_param(func, init): 
	if init:
		paramstr=check_output(['bin/read_params', 'init-param/'+func])
	else:
		paramstr=check_output(['bin/read_params', 'param/'+func])
	paramv=paramstr.split(' ')
	paramv.pop()
	params = []
	for elem in paramv:	
		params += [float(elem)]
	return params

#	mult_energy_ - Get array of energies from Gauss	
#########################################################
#	Input:
#	systems - Name of the systems
#########################################################
#	Output:
#	E_mult - dictionary of energies
#########################################################
def mult_energy(systems):
	E_mult = {}
	for system in systems:
		E_mult[system] = [get_energy(system)]
	return E_mult

#	run_mult - Run a set of jobs 
#########################################################
#	Input:
#	system - The name of the system
#########################################################
def run_mult(systems):
	for system in systems:
		os.system('./run ' + system)

#	wait_mult - Wait until a set of jobs are done	
#########################################################
#	Input:
#	systems - The names of the systems
#########################################################
def wait_mult(systems):
	for system in systems:
		while bool(int(check_output(['bin/number_of_runs', system]))):
			pass

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
#	Input:
#	array - Array to be saved
#	filename - Path to and name of file
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
#	Input:
#	filename - Path to and name of file
#########################################################
#	Output:
#	array - Array that's loaded
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

#	write - write a string to the output file
#########################################################
#	Input:
#	instr - String to be written
#########################################################
def write(instr):
	os.system("echo '"+instr+"' >> output")

#	read_bench - Reads benchmark energy from .ener file
#########################################################
#	Input:
#	system - system to get benchmark
#########################################################
#	Output:
#	bench - A dictionary containing info for this
#	particular system
#########################################################
def read_bench(system):
	bench_data = load('data/ener/'+system+'.ener')
	bench = {
		# Intermolecular distance
		#'R': bench_data[0], 
		# Total energy
		'E': bench_data[-10],
		# Monomer A energy
		#'E_A': bench_data[-2],
		# Monomer B energy
		#'E_B': bench_data[-1]
	}
	return bench

#	read_bench_mult 
#########################################################
#	Input:
#	systems - array of systems to read benchmark 
#	data from
#########################################################
#	Output:
#	bench - dictionary of all benchmark data
#########################################################
def read_bench_mult(systems):
	bench = {}
	for system in systems:
		bench[system] = read_bench(system)
	return bench

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
