from subprocess import check_output
import os, pdb, re

# data structure containing all the important information 
# such as energy for each system from each file
# systems are referred to by by NAMEKEY
# datatype of that system is referred to as TYPEKEY
class data:

	# NAMEKEYS - list of NAMEKEYS used
	def __init__(self, NAMEKEYS=[]):
		self.iteration = 0
		self.NAMEKEYS = NAMEKEYS	
		self.E = {}
		self.error = 0
	
		for NAMEKEY in NAMEKEYS:
			self.E[NAMEKEY] = {}
				
	# run - run all NAMEKEYS of specified input file
	def run(self, TYPEKEY): 
		for NAMEKEY in self.NAMEKEYS:
			os.system('bin/run '+TYPEKEY+' '+NAMEKEY)

	# wait - wait for all NAMEKEYS being ran
	def wait(self):
		for NAMEKEY in self.NAMEKEYS:
			while bool(int(check_output(['bin/number_of_runs', NAMEKEY]))):
                        	pass

	# energy - get energies from specified file type
	def energy(self, TYPEKEY):
		for NAMEKEY in self.NAMEKEYS:
			try:
				self.E[NAMEKEY][TYPEKEY]=float(check_output(['bin/energy_from_'+TYPEKEY, NAMEKEY]))
			except:
				print('Could not get '+TYPEKEY+' energy for '+NAMEKEY)
				self.E[NAMEKEY][TYPEKEY]=False		

	# convert - converts file TYPEKEY_IN to file TYPEKEY_OUT
	def convert(self, TYPEKEY_IN, TYPEKEY_OUT):
		for NAMEKEY in self.NAMEKEYS:
			os.system('bin/'+TYPEKEY_IN+'2'+TYPEKEY_OUT+' '+NAMEKEY)

	#mult - multiplies all data of a certain type by a constant x
	# useful for unit conversion and flipping signs
	def mult(self, TYPEKEY, x):
		for NAMEKEY in self.NAMEKEYS:
			self.E[NAMEKEY][TYPEKEY] = x*self.E[NAMEKEY][TYPEKEY]

	# update_error - add up all energies in TYPEKEYS, and divide by DENOM. 
	# be sure to change the signs beforehand by using the times function
	def update_error(self, DENOM, *TYPEKEYS):
		self.error=0
		for NAMEKEY in self.NAMEKEYS:
			self.error += sum([E[NAMEKEY][TYPEKEY] for TYPEKEY in TYPEKEYS])/E[NAMEKEY][DENOM]

	# read_input - reads NAMEKEYS from input file
	def read_input(self):
		f = open('input', 'r')
		f.seek(0)
		fstr = f.read()
		fvect = fstr.split('\n')
		fvect.pop()
		f.close()
		for NAMEKEY in fvect:
			self.E[NAMEKEY] = {}
			self.NAMEKEYS += [NAMEKEY]
	
	# write_output - a preset message to the output file, giving dates and energies, also calculates error
	def write_output(self):
		w = lambda s : os.system("echo '"+s+"' >> output")
		w('\n')
		os.system('echo $(date + "%H:%M) >> output')
		w('Iteration: '+str(self.iteration))
		for NAMEKEY in self.NAMEKEYS:
			for TYPEKEY in self.E[NAMEKEY]:
				try:
					w(TYPEKEY+' energy for '+NAMEKEY+' : '+E[NAMEKEY][TYPEKEY])
				except: 
					print('Failed to load '+NAMEKEY+'.'+TYPEKEY+' energy')
		w('-------------------------')
		w('Total error : '+self.error)
		w('-------------------------')	

# save_array - save an array to a filename, space seperated
def save_array(array, filename):
        fstr = ''
        for i in range(len(array)):
                fstr += str(array[i]) + ' '
        fstr += '\n'
        f = open(filename, 'w')
        f.seek(0)
        f.write(fstr)
        f.close()

# load_array - load an array from a filename, space seperated
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
