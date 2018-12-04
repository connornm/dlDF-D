from subprocess import check_output
import os, pdb, re

class data:


	# systems - list of systems used
	# test - filetype for software being used for the test (com fo Gaussion)
	# bench - filetype for software being used as the benchmark (ener for autoPES) 
	# disp - filetype for software used to calculate just dispersion energy (cnf for disp4dldf.py)
	# loadin - boolean to signal autmoatic load systems from a file named 'input' in current directory
	def __init__(self, systems=[], loadin=False):
		self.iteration = 0
		self.systems = systems	
		self.E = {}
	
		for system in systems:
			self.E[system] = {}
	
		if loadin:
			f = open('input', 'r')
			f.seek(0)
			fstr = f.read()
			fvect = fstr.split('\n')
			fvect.pop()
			f.close()
			for elem in fvect:
				self.systems += [elem]


	# write - write to the output
	def write(self):
		w = lambda s : os.system("echo '"+s+"' >> output")
		w('\n')
		os.system('echo $(date + "%H:%M) >> output')
		w('Iteration: '+str(self.iteration))
		for system in self.systems:
			for datatype in self.E[system]:
				try:
					w(datatype+' energy for '+system+' : '+E[system][datatype])
				except: 
					print('Failed to load '+system+'.'+datatype+' energy')
			
	# run - run all systems of specified input file
	def run(self, datatype): 
		for system in self.systems:
			os.system('./run'+datatype+' '+system)

	# wait - wait for all systems being ran
	def wait(self):
		for system in self.systems:
			while bool(int(check_output(['bin/number_of_runs', system]))):
                        	pass

	# energy - get energies from specified file type
	def energy(self, datatype):
		for system in self.systems:
			try:
				self.E[system][datatype]=float(check_output(['bin/energy_from_'+datatype, self.system]))
			except:
				print('Could not get '+datatype+' energy for '+self.system)
				self.E[system][datatype]=False			

		
	# xyz2com - convert all xyz files to com files
	def xyz2com(self):
		for system in self.systems:
			os.system('bin/xyz2com '+system)

	# com2cnf - convert all com files to cnf files
	def com2cnf(self):
		for system in self.systems:
			os.system('bin/com2cnf '+system)


