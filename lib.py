from subprocess import check_output
import os, pdb, re, pickle

# data structure containing all the important information 
# such as energy for each system from each file
# systems are referred to by by NAMEKEY
# datatype of that system is referred to as TYPEKEY
class data:

	# NAMEKEYS - list of NAMEKEYS used
	def __init__(self, QUANKEYS=[], NAMEKEYS=[], TYPEKEYS=[]):
		self.Hartree = 627.5095 #kcal/mol
		self.kcalpmol = 1/self.Hartree #Hartrees
		self.iteration = 0
		self.error = 0	
		self.NAMEKEYS = NAMEKEYS
		self.QUANKEYS = QUANKEYS
		self.TYPEKEYS = TYPEKEYS
		self.vals = {}

	# run - run all NAMEKEYS of specified TYPEKEY input file
	def run(self, TYPEKEY): 
		for NAMEKEY in self.NAMEKEYS:
			os.system('bin/run '+TYPEKEY+' '+NAMEKEY)

	# wait - wait for all NAMEKEYS being ran
	def wait(self):
		for NAMEKEY in self.NAMEKEYS:
			while bool(int(check_output(['bin/number_of_runs', NAMEKEY]))):
                        	pass

	# get - get all of quantity from certain type
	def get(self, QUANKEY, TYPEKEY):
		if TYPEKEY not in self.TYPEKEYS:
			self.TYPEKEYS += [TYPEKEY]
		if QUANKEY not in self.QUANKEYS:
			self.QUANKEYS += [QUANKEY]
		for NAMEKEY in self.NAMEKEYS:
			try:
				self.vals[tuple((NAMEKEY,TYPEKEY,QUANKEY))]=float(check_output(['bin/'+QUANKEY+'_from_'+TYPEKEY, NAMEKEY]))
			except:
				print('Could not get '+TYPEKEY+' '+QUANKEY+' for '+NAMEKEY)
				self.vals[tuple((NAMEKEY,TYPEKEY,QUANKEY))]=None	

        #multiplys - multiplies all data of a certain type by a constant x
        # useful for unit conversion and flipping signs
	def multiply(self, TYPEKEY, QUANKEY, x):
		for NAMEKEY in self.NAMEKEYS:
			self.vals[tuple((NAMEKEY,TYPEKEY,QUANKEY))] = x*self.vals[(NAMEKEY, TYPEKEY, QUANKEY)]

	# add - makes a new type by adding other energies together
	def add(self, NEWTYPEKEY, QUANKEY, *TYPEKEYS):
		for NAMEKEY in self.NAMEKEYS:
			self.vals[tuple((NAMEKEY,NEWTYPEKEY,QUANKEY))]=sum([self.vals[(NAMEKEY,TYPE,QUANKEY)] for TYPE in TYPEKEYS])

	# convert - converts file TYPEKEY_IN to file TYPEKEY_OUT
	def convert(self, TYPEKEY_IN, TYPEKEY_OUT):
		for NAMEKEY in self.NAMEKEYS:
			if os.system('bin/'+TYPEKEY_IN+'2'+TYPEKEY_OUT+' '+NAMEKEY):
				return KeyError	
			
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
			self.NAMEKEYS += [NAMEKEY]
	
	# write_output - writes the energies in types in order given
	def write_output(self, DTYPE, *TYPEKEYS):
		def wr(s):
			os.system("echo '"+s+"' >> output")
		for NAMEKEY in self.NAMEKEYS:
			out=str(self.vals[(NAMEKEY,DTYPE,'distance')])
			for i in range(len(TYPEKEYS)):
				out += ' '+str(self.vals[(NAMEKEY,TYPEKEYS[i],'energy')])	
			wr(out)

	# save - saves current self to f
	def save(self, f='quick'):
		f='.save/'+f
		pickle.dump(self.vals, open(f+'vals', 'wb'))
		pickle.dump(self.NAMEKEYS, open(f+'name', 'wb'))
		pickle.dump(self.TYPEKEYS, open(f+'type', 'wb'))
		pickle.dump(self.QUANKEYS, open(f+'quan', 'wb'))

	# load - loads self from f
	def load(self, f='quick'):
		f='.save/'+f
		self.vals=pickle.load(open(f+'vals', 'rb'))
		self.NAMEKEYS=pickle.load(open(f+'name', 'rb'))
		self.TYPEKEYS=pickle.load(open(f+'type', 'rb'))
		self.QUANKEYS=pickle.load(open(f+'quan', 'rb'))


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



