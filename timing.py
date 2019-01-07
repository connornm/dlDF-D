import sys, random, time
sys.path.append('praxis/')
from praxis import praxis
import lib, os, pdb
tolerance = 0.03
max_step = 1
from numpy import log, exp
software='gaussian'

# Write to the output
def w(s):
	os.system("echo '"+s+"' >> output")	

data = lib.data()
data.meta['system']={('qqq'+str(2*i).zfill(5)) for i in range(1, 15)}
data.convert('xyz', software)
data.run(software)
data.wait('qqq')
data.get('time', software)
atoms, times = [], []
for system in data.meta['system']:
	atoms += [int(system[4:])]
	times += [data.vals[system][software]['time']]

print(atoms)
print(times)
