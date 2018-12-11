from subprocess import check_output
import os, pdb, re

# data structure containing all the important information 
# such as energy for each system from each file
# systems are referred to by by NAMEKEY
# datatype of that system is referred to as TYPEKEY
class data:

	# NAMEKEYS - list of NAMEKEYS used
	def __init__(self, NAMEKEYS=[]):
		self.Hartree = 627.5095 #kcal/mol
		self.kcalpmol = 1/self.Hartree #Hartrees
		self.iteration = 0
		self.error = 0	
		self.NAMEKEYS = NAMEKEYS	
		self.E = {}
		self.R = {}	
		for NAMEKEY in NAMEKEYS:
			self.E[NAMEKEY] = {}
			self.R[NAMEKEY] = {}	

	# run - run all NAMEKEYS of specified TYPEKEY input file
	def run(self, TYPEKEY): 
		for NAMEKEY in self.NAMEKEYS:
			os.system('bin/run '+TYPEKEY+' '+NAMEKEY)

	# wait - wait for all NAMEKEYS being ran
	def wait(self):
		for NAMEKEY in self.NAMEKEYS:
			while bool(int(check_output(['bin/number_of_runs', NAMEKEY]))):
                        	pass

	# get_distance - get COM-COM distances
	def get_distance(self, TYPEKEY):
		for NAMEKEY in self.NAMEKEYS:
			try:
				self.R[NAMEKEY][TYPEKEY]=float(check_output(['bin/distance_from_'+TYPEKEY, NAMEKEY]))
			except:
				print('Could not get '+TYPEKEY+' distance for '+NAMEKEY)
				self.R[NAMEKEY][TYPEKEY]=False		


	# get_energy - get energies from specified file type
	def get_energy(self, TYPEKEY):
		for NAMEKEY in self.NAMEKEYS:
			try:
				self.E[NAMEKEY][TYPEKEY]=float(check_output(['bin/energy_from_'+TYPEKEY, NAMEKEY]))
			except:
				print('Could not get '+TYPEKEY+' energy for '+NAMEKEY)
				self.E[NAMEKEY][TYPEKEY]=False		

        #multiplys - multiplies all data of a certain type by a constant x
        # useful for unit conversion and flipping signs
	def multiply(self, TYPEKEY, x):
		for NAMEKEY in self.NAMEKEYS:
			self.E[NAMEKEY][TYPEKEY] = x*self.E[NAMEKEY][TYPEKEY]

	# add - makes a new type by adding other energies together
	def add(self, NEWTYPEKEY, *TYPEKEYS):
		for NAMEKEY in self.NAMEKEYS:
			self.E[NAMEKEY][NEWTYPEKEY]=sum([self.E[NAMEKEY][TYPE] for TYPE in TYPEKEYS])

	# convert - converts file TYPEKEY_IN to file TYPEKEY_OUT
	def convert(self, TYPEKEY_IN, TYPEKEY_OUT):
		for NAMEKEY in self.NAMEKEYS:
			os.system('bin/'+TYPEKEY_IN+'2'+TYPEKEY_OUT+' '+NAMEKEY)
	
	# read_input - reads NAMEKEYS from input file
	def read_input(self):
		self.NAMEKEYS = []
		f = open('input', 'r')
		f.seek(0)
		fstr = f.read()
		fvect = fstr.split('\n')
		fvect.pop()
		f.close()
		for NAMEKEY in fvect:
			self.E[NAMEKEY] = {}
			self.R[NAMEKEY] = {}
			self.NAMEKEYS += [NAMEKEY]
	
	# write_output - writes the energies in types in order given
	def write_output(self, DTYPE, *TYPEKEYS):
		def wr(s):
			os.system("echo '"+s+"' >> output")
		for NAMEKEY in self.NAMEKEYS:
			out=str(self.R[NAMEKEY][DTYPE])
			for i in range(len(TYPEKEYS)):
				out += ' '+str(self.E[NAMEKEY][TYPEKEYS[i]])	
			wr(out)

# new_list - make a new list of inputs
def new_list(SYSTEM, NUMBER):
	os.system('bin/new_list '+SYSTEM+' '+str(NUMBER))

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
def load_array(filename):
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



