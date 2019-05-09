#Various useful functions stored here to reduce clutter in the main loop
import sys, random, os, pdb, pickle
from numpy import log, exp, sqrt


# Write to the output
def w(s, raw=False):
	if raw:
		os.system("echo '"+s+"' >> raw_output")
	else:
		os.system("echo '"+s+"' >> output")	

# Gets the error
def get_error(dat):
	error = 0
	E_mon = dat.vals['BLIND_monomer']['gaussian']['energy']	
#	E_mon = dat.vals['BLIND_monomer']['gaussian']['energy']	
	os.system('rm -f raw_output')
	for system in dat.meta['dimer']:	
#		E_monA = dat.vals[system+'_mA']['gaussian']['energy']
#		E_monB = dat.vals[system+'_mB']['gaussian']['energy']
		E_dim = dat.vals[system]['gaussian']['energy']
		E_disp = dat.vals[system]['disp']['energy']
		E_sapt = dat.vals[system]['sapt']['energy']
		w(system)
		if E_dim and E_mon and E_disp and E_sapt:
#			dE = dat.Hartree*(E_dim-E_monA-E_monB)+E_disp
			dE = dat.Hartree*(E_dim-2*E_mon) + E_disp
			er = (E_sapt-dE)**2
			w('All energies in kcal/mol')
			w('E_dim='+str(E_dim))
#			w('E_monA='+str(E_monA))
#			w('E_monB='+str(E_monB))
			w('E_mon='+str(E_mon))
			w('E_disp='+str(E_disp))		
			w('E_sapt='+str(E_sapt))	
			w('E_dim-2*E_mon+E_disp='+str(dE))
			w('er='+str(er))
			w('"'+system+'"'+':{"E_sapt":'+str(E_sapt)+', "er":'+str(er)+'}, ', raw=True)
			error += er
		else:
			w(system)
			w('Failed to converge')
			w('Adding 1000000 to error')
			error += 1000000
		w(' ')
	return error

def print_param(param, error, iteration):
	w(' ')
	os.system('echo $(date) >> output')
	w('Iteration '+str(iteration))
	w('Total Error: '+str(error))
	w('Parameters:')
	for i in range(len(param)):
		w(str(param[i])+',')
	w('--------------------------------------------')
	w(' ')

# save_array
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

# sigmoid function:
def sig(x):
	return 1.0/(exp(x)+1.0)

# inverse sigmoid function
def isig(x):
	return log(1.0/x - 1.0 )
