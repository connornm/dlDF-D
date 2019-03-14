from subprocess import check_output
import os, pdb, re, pickle

# auto-filling dictionary
class vdict(dict):
	def __init__(self):
		self={}
	def __missing__(self,key):
		self[key]=vdict()
		return self[key]

# data structure containing all the important information 
class data:

	# NAMEKEYS - list of NAMEKEYS used
	def __init__(self):
		self.Hartree = 627.5095 #kcal/mol
		self.kcalpmol = 1/self.Hartree #Hartrees
		self.iteration=0
		# holds all the values for data 
		self.vals=vdict()
		# metadata for names of systems, softwares, and quantities used
		self.meta={
			'system': set(),
			'software': set(),
			'quantity': set(),
			'user': {'connornm'}
			}

	# split - splits all molecules in to two molecules at line number
	# supports user defined function f for naming convention
	def split(self, software_in, software_out, n, system_in, system_out1, system_out2):
		self.meta['software'].add(software_in), self.meta['system'].add(system_out1), self.meta['system'].add(system_out2)
		command='bin/split_'+software_in+'2'+software_out+'_'+software_out
		os.system(command+' '+system_in+' '+str(n)+' '+system_out1+' '+system_out2)
		
	# merge - adds two molecules into a single molecule
	def merge(self, software_in, software_out, system_in1, system_in2, system_out):
		self.meta['software'].add(software_in), self.meta['system'].add(system_out)
		command='bin/merge_'+software_in+'_'+software_in+'2'+software_out
		os.system(command+' '+system_in1+' '+system_in2+' '+system_out)

	# run - run all the systems for a certain software
	def run(self, software, key='system'): 
		self.meta['software'].add(software)
		for system in self.meta[key]:
			os.system('bin/run '+software+' '+system)

	# wait - wait for all NAMEKEYS being ran
	def wait(self, string):
		while bool(int(check_output(['bin/number_of_runs', string]))):
                        pass

	# get - get all of quantity from certain software
	def get(self, quantity, software, key='system'):	
		self.meta['software'].add(software)
		self.meta['quantity'].add(quantity)
		for system in self.meta[key]:
			try:
				raw=check_output(['bin/'+quantity+'_from_'+software, system])
				self.vals[system][software][quantity]=float(raw)
			except:
				print('Could not get '+software+' '+quantity+' for '+system)
				self.vals[system][software][quantity]=None	

	# convert - converts file software_in to file software_out
	def convert(self, software_in, software_out, key='system'):
		for system in self.meta[key]:
			if os.system('bin/'+software_in+'2'+software_out+' '+system):
				return KeyError	

	# save - saves current self to file
	def save(self, f='quick'):
		f='.save/'+f
		pickle.dump(self.vals, open(f+'vals', 'wb'))
		pickle.dump(self.meta, open(f+'meta', 'wb'))

	# load - loads self from f
	def load(self, f='quick'):
		f='.save/'+f
		self.vals=pickle.load(open(f+'vals', 'rb'))
		self.meta=pickle.load(open(f+'meta', 'rb'))


