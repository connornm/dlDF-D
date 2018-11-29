#!/usr/bin/env python
#
# This Python script calculates the dispersion part of the interaction 
# energy between two systems using the expression fitted in Podeszwa
# et al., JPCL 1, 550 (2010) (version=1) or the expression fitted in
# November 2010, to be published (version=2). The version is specified
# as the second parameter in the command line.
#
# Both numerical and analytical gradients of dispersion with respect to
# nuclear coordinates are available.
#
# Usage:
#
# 1. disp4dldf.py input_file version optional_arguments
# Read the geometry from input_file.
# input_file may be a Gaussian input file (with some restrictions) or
# can roughly follow the dimer.cnf file format from SAPT:
#
# empty line(s)
# X(atom 1) Y(atom 1) Z(atom 1) charge(atom 1) .... \
# X(atom 2) Y(atom 2) Z(atom 2) charge(atom 2) ....  \ monomer
# X(atom 3) Y(atom 3) Z(atom 3) charge(atom 3) ....  /    A
#   ...       ...       ...          ...       .... /
# empty line(s)      
# X(atom 1) Y(atom 1) Z(atom 1) charge(atom 1) .... \
# X(atom 2) Y(atom 2) Z(atom 2) charge(atom 2) ....  \ monomer
# X(atom 3) Y(atom 3) Z(atom 3) charge(atom 3) ....  /    B
#   ...       ...       ...          ...       .... /
# empty line(s)      
# whatever
#
# where all the coordinates X,Y,Z are given in atomic units.
# The 'empty line(s)' need not be completely empty, but neither of them
# should contain more than three fields separated by spaces.
#
# In this mode, the following optional arguments can be supplied, in any
# order after the first two: 
# 1 - a single subsystem is assumed and the dispersion energy of the whole
# system is calculated.
# n - compute numerical gradients with respect to nuclear coordinates.
# a - compute analytical gradients with respect to nuclear coordinates.
#
# 2. disp4dldf.py dir_with_saptdft_outputs version
# Read all .out files from the directory dir_with_saptdft_outputs, treating 
# them as SAPT(DFT) output files. For each file, extract geometry and 
# the benchmark result, calculate the fitted dispersion energy, and make
# a nice comparison with MUE and MURE.
#
# At short distances, dispersion is switched off using the Fermi function
# 1.0/(1.0+math.exp(-bswitch*(r/r0-1.0)))
# where
# r0=aswitch*(covalent_radius(at1)+covalent_radius(at2))
# and aswitch,bswitch are parameters read from the switch.par file in the
# current directory (if the file is not there, reasonable default values
# are chosen).
# 
import math
import os
import sys

def catomatom(at1,at2,dict):
  """Returns C6 coefficient between atoms with charges at1 and at2."""
  c6a=catom(at1,dict)
  c6b=catom(at2,dict)
  return math.sqrt(c6a*c6b)

def catom(at,dict):
  """Returns C6 coefficient for atom at, converted to kcal*bohr^6/mol."""
  c=dict[at]
  return c*10884.3258684610
# the coefficient is 1/(4184.0*0.0529177209^6)

def catomatom8(at1,at2,dict8):
  """Returns C8 coefficient between atoms with charges at1 and at2."""
  c8a=catom8(at1,dict8)
  c8b=catom8(at2,dict8)
  return math.sqrt(c8a*c8b)

def catom8(at,dict8):
  """Returns C8 coefficient for atom at, converted to kcal*bohr^8/mol."""
  c=dict8[at]
  return c*3886863.35457252
# the coefficient is 1/(4184.0*0.0529177209^8)

def damp_TT(r,at1,at2,dicd):
  """Returns the value of Tang-Toennies damping function at separation r."""
# it is assumed below that the damping factor is the geometric mean of atomic
# damping factors:
  alpha=math.sqrt(abs(dicd[at1]*dicd[at2]))
  br=alpha*r
  sum=1.0
  term=1.0
  for i in [1,2,3,4,5,6]:
    term=term*br/i
    sum=sum+term
  d=1.0-math.exp(-br)*sum
  return d

def damp_TT8(r,at1,at2,dicd):
  """Returns the value of Tang-Toennies damping function f_8 at separation r."""
# it is assumed below that the damping factor is the geometric mean of atomic
# damping factors:
  alpha=math.sqrt(abs(dicd[at1]*dicd[at2]))
  br=alpha*r
  sum=1.0
  term=1.0
  for i in [1,2,3,4,5,6,7,8]:
    term=term*br/i
    sum=sum+term
  d=1.0-math.exp(-br)*sum
  return d

def covalent_radius(at):
  """Returns the covalent (Bragg) radius of the atom at, converted to bohr.
  The radii are taken from the Cambridge Structural Database."""
  try:
    dummy=len(at)
  except:
    if (at==2):
      r=1.50  
    if (at==3):
      r=1.28  
    if (at==4):
      r=0.96  
    if (at==5):
      r=0.83  
    if (at==6):
      r=0.68  
    if (at==7):
      r=0.68  
    if (at==8):
      r=0.68  
    if (at==9):
      r=0.64  
    if (at==10):
      r=1.50  
    if (at==11):
      r=1.66  
    if (at==12):
      r=1.41  
    if (at==13):
      r=1.21  
    if (at==14):
      r=1.20  
    if (at==15):
      r=1.05  
    if (at==16):
      r=1.02  
    if (at==17):
      r=0.99  
    if (at==18):
      r=1.51  
  else:
    if at[0]==1:
      r=0.23
    elif at[0]==3:
      r=1.28  
    elif at[0]==4:
      r=0.96  
    elif at[0]==11:
      r=1.66  
    elif at[0]==12:
      r=1.41  
    else:
      print ("Wrong atom type in covalent_radius")
      raise "Error!"
  return r/0.529177209  

def switching_function(r,at1,at2,aswitch,bswitch):
  """Returns the value of the dispersion switching (Fermi) function for
  atoms at1,at2 at separation r."""
  r0=aswitch*(covalent_radius(at1)+covalent_radius(at2))
  swf=1.0/(1.0+math.exp(-bswitch*(r/r0-1.0)))
  return swf

def disp_total(r,at1,at2,dict,dict8,dicd,aswitch,bswitch):
  """Returns the overall contribution to dispersion 
  energy (in kcal/mol) for atoms at1,at2 at separation r."""
  if not(at1 in dict.keys()):
    print ("Dispersion coefficients not defined for atom",at1)
    raise "Error!"
  if not(at2 in dict.keys()):
    print ("Dispersion coefficients not defined for atom",at2)
    raise "Error!"
  edisp=-catomatom(at1,at2,dict)/(r*r*r*r*r*r)*damp_TT(r,at1,at2,dicd)
  edisp-=catomatom8(at1,at2,dict8)/(r*r*r*r*r*r*r*r)*damp_TT8(r,at1,at2,dicd)
  edisp*=switching_function(r,at1,at2,aswitch,bswitch)
  return edisp

def read_saptdft_output(loutput,version):
  """Extract geometry information and benchmark (disp+exch-disp) result from
  the SAPT(DFT) output that is supplied as a list of lines."""
  geodone=0
  geocheck=0
  imode=0
  icount=0
  eexdispcks=0.0
  laa=[]
  lab=[]
  for s1 in loutput:
    s1=s1.strip("\n")
    ll=s1.split()
    if (len(ll)==0):
      imode=0
      continue
    if (ll[0]=="-----------"):
      if (len(ll)<3):
        continue
      if (ll[1]=="Molecule") and (ll[2]=="A"):
        imode=1
        icount=0
        if (geodone>=1):
          geocheck=1
          laaold=laa
          labold=lab
          laa=[]
          lab=[]
      if (ll[1]=="Molecule") and (ll[2]=="B"):
        imode=2
        geodone=geodone+1
        icount=0
      continue
    if (ll[0]=="ATOM"):
      continue
    if (imode==1):
      atom=(float(ll[1]),float(ll[2]),float(ll[3]),int(float(ll[4])))
      laa.append(atom)
      if (geocheck==1) and (laa[icount]!=laaold[icount]):
        print ("Inconsistent geometry for A, count=",icount,"\n")
      icount=icount+1
      continue
    if (imode==2):
      if (len(ll)<5):
        imode=0
        continue
      atom=(float(ll[1]),float(ll[2]),float(ll[3]),int(float(ll[4])))
      lab.append(atom)
      if (geocheck==1) and (lab[icount]!=labold[icount]):
        print ("Inconsistent geometry for B, count=",icount,"\n")
      icount=icount+1
      continue
    if (ll[0]=="E^{(20)}_{disp}"):
      edispucks=float(ll[2])
      continue
    if (ll[0]=="E^{(20)}_{exch-disp}"):
      eexdispucks=float(ll[2])
      continue
    if (len(ll)<3):
      continue
    if (ll[0]=="CKS") and (ll[1]=="dispersion") and (ll[2]!="and"):
      edispcks=float(ll[3])
    if (ll[0]=="CKS") and (ll[1]=="dispersion:") and (ll[2]!="and"):
      edispcks=float(ll[3])
    if (ll[0]=="Final") and (ll[1]=="CKS") and (ll[2]=="dispersion"):
      edispcks=float(ll[5])
    if (len(ll)<6):
      continue
    if ll[0]=="e200" and ll[1]=="exch" and ll[3]=="disp" and ll[4]=="(CKS)":
      if ll[5]==":":
        eexdispcks=float(ll[6].replace("D","E"))*627.51
      else:
        eexdispcks=float((ll[5][1:]).replace("D","E"))*627.51
  if geodone==0:
# no geometry in summary table for open-shell SAPT(DFT). We have to extract
# the geometry from the DALTON output.
    geostart=0
    geoend=0
    for s1 in loutput:
      s1=s1.strip("\n")
      ll=s1.split()
      if geostart==0:
        if len(ll)<4:
          continue
        if ll[0]=="Total" and ll[1]=="number" and ll[3]=="coordinates:":
          geostart=1
      elif geoend==1:
        if len(ll)<4:
          continue
        if ll[0]=="Dispersion(CKS):":
          edispcks=float(ll[2])
        if ll[0]=="Exchange-dispersion(CKS)*:":
          eexdispcks=float(ll[2])
      else:
        if len(ll)<3:
          continue
        if ll[0]=="Interatomic":
          geoend=1
          continue
        if ll[1]=="Mb":
          geostart=2
        elif ll[2]=="x":
          if geostart==2:
            geostart=3
          x=float(ll[3])
          nat=find_atom_type(ll[1])
        elif ll[1]=="y":
          y=float(ll[2])
        elif ll[1]=="z":  
          z=float(ll[2])
          atom=(x,y,z,nat)
          if geostart==1:
            laa.append(atom)
          elif geostart==3:
            lab.append(atom)
  if eexdispcks==0.0:
    eexdispcks=eexdispucks*(edispcks/edispucks)
  edisp=edispcks
  eexdisp=eexdispcks
  res=edisp+eexdisp
  ltemp=laa
  laa=analyze_hydrogens(ltemp)
  ltemp=lab
  lab=analyze_hydrogens(ltemp)
  if version==2:
    ltemp=laa
    laa=analyze_atoms2(ltemp,3)
    ltemp=lab
    lab=analyze_atoms2(ltemp,3)
    ltemp=laa
    laa=analyze_atoms(ltemp,4)
    ltemp=lab
    lab=analyze_atoms(ltemp,4)
    ltemp=laa
    laa=analyze_atoms2(ltemp,11)
    ltemp=lab
    lab=analyze_atoms2(ltemp,11)
    ltemp=laa
    laa=analyze_atoms(ltemp,12)
    ltemp=lab
    lab=analyze_atoms(ltemp,12)
  return (laa,lab,res)

def read_geometry(li,version):
  """Extract geometry information from the list of lines (dimer.cnf format)."""
  imode=0
  isep=0
  laa=[]
  lab=[]
  for s in li:
    if (imode>2):
      break
    s=s.strip("\n")
    ll=s.split()
    if (len(ll)<4):
      if (isep==0):
        imode=imode+1
      isep=1
      continue
    isep=0
    atom=(float(ll[0]),float(ll[1]),float(ll[2]),int(float(ll[3])))
    if (imode==1):    
      laa.append(atom)
    else:
      lab.append(atom)
  ltemp=laa
  laa=analyze_hydrogens(ltemp)
  ltemp=lab
  lab=analyze_hydrogens(ltemp)
  if version==2:
    ltemp=laa
    laa=analyze_atoms2(ltemp,3)
    ltemp=lab
    lab=analyze_atoms2(ltemp,3)
    ltemp=laa
    laa=analyze_atoms(ltemp,4)
    ltemp=lab
    lab=analyze_atoms(ltemp,4)
    ltemp=laa
    laa=analyze_atoms2(ltemp,11)
    ltemp=lab
    lab=analyze_atoms2(ltemp,11)
    ltemp=laa
    laa=analyze_atoms(ltemp,12)
    ltemp=lab
    lab=analyze_atoms(ltemp,12)
  return (laa,lab)

def analyze_hydrogens(latom):
  """Find all hydrogens in the molecule and label them according to their
  nearest neighbor."""
  latomout=latom
  for i in range(len(latom)):
    atom=latom[i]
    atomout=latomout[i]
    if (atom[3]!=1):
      continue
    closestr=100.0
    for atom2 in latom:
      if (atom2==atom):
        continue
      r=math.sqrt((atom[0]-atom2[0])*(atom[0]-atom2[0])+(atom[1]-atom2[1])*(atom[1]-atom2[1])+(atom[2]-atom2[2])*(atom[2]-atom2[2]))
      if (r<closestr):
        closestr=r
        type=atom2[3]
    try:
      type=type[0]
    except:
      pass
    if (type==1):
#     print "H connected to H"
      atomout=(atom[0],atom[1],atom[2],(1,1))
    if (type==3):
#     print "H connected to Li"
      atomout=(atom[0],atom[1],atom[2],(1,3))
    if (type==4):
#     print "H connected to Be"
      atomout=(atom[0],atom[1],atom[2],(1,4))
    if (type==5):
#     print "H connected to B"
      atomout=(atom[0],atom[1],atom[2],(1,5))
    if (type==6):
#     print "H connected to C"
      atomout=(atom[0],atom[1],atom[2],(1,6))
    if (type==7):
#     print "H connected to N"
      atomout=(atom[0],atom[1],atom[2],(1,7))
    if (type==8):
#     print "H connected to O"
      atomout=(atom[0],atom[1],atom[2],(1,8))
    if (type==9):
#     print "H connected to F"
      atomout=(atom[0],atom[1],atom[2],(1,9))
    if (type==11):
#     print "H connected to Na"
      atomout=(atom[0],atom[1],atom[2],(1,11))
    if (type==12):
#     print "H connected to Mg"
      atomout=(atom[0],atom[1],atom[2],(1,12))
    if (type==13):
#     print "H connected to Al"
      atomout=(atom[0],atom[1],atom[2],(1,13))
    if (type==14):
#     print "H connected to Si"
      atomout=(atom[0],atom[1],atom[2],(1,14))
    if (type==15):
#     print "H connected to P"
      atomout=(atom[0],atom[1],atom[2],(1,15))
    if (type==16):
#     print "H connected to S"
      atomout=(atom[0],atom[1],atom[2],(1,16))
    if (type==17):
#     print "H connected to Cl"
      atomout=(atom[0],atom[1],atom[2],(1,17))
    if not(type in [1,3,4,5,6,7,8,9,11,12,13,14,15,16,17]):
      print ("Wrong neighbor! ",type)
      raise "Error"
    latomout[i]=atomout
  return latomout

def analyze_atoms(latom,natom):
  """Find all atoms "natom" in the molecule and label them as
     (natom,0) - no neighbors (natom,1) - one neighbor ..."""
# use this for Be,Mg
  latomout=latom
  for i in range(len(latom)):
    atom=latom[i]
    atomout=latomout[i]
    if (atom[3]!=natom):
      continue
    nneighbors=0
    for atom2 in latom:
      if (atom2==atom):
        continue
      r=math.sqrt((atom[0]-atom2[0])*(atom[0]-atom2[0])+(atom[1]-atom2[1])*(atom[1]-atom2[1])+(atom[2]-atom2[2])*(atom[2]-atom2[2]))
      if (r<3.9): # constant OK for Li,Be,Na,Mg
        nneighbors=nneighbors+1
    atomout=(atom[0],atom[1],atom[2],(natom,nneighbors))
    latomout[i]=atomout
  return latomout

def analyze_atoms2(latom,natom):
  """Find all atoms "natom" in the molecule and label them as
     (natom,0) - no neighbors (natom,1) - neighbor H
     (natom,2) - other neighbor."""
# use this for Li,Na
  latomout=latom
  for i in range(len(latom)):
    atom=latom[i]
    atomout=latomout[i]
    if (atom[3]!=natom):
      continue
    nneighbors=0
    for atom2 in latom:
      if (atom2==atom):
        continue
      r=math.sqrt((atom[0]-atom2[0])*(atom[0]-atom2[0])+(atom[1]-atom2[1])*(atom[1]-atom2[1])+(atom[2]-atom2[2])*(atom[2]-atom2[2]))
      if (r<3.9): # constant OK for Li,Be,Na,Mg
        nneighbors=nneighbors+1
        type=atom2[3]
    try:
      type=type[0]
    except:
      pass
    if nneighbors==0:
      atomout=(atom[0],atom[1],atom[2],(natom,0))
    elif (type==1):
      atomout=(atom[0],atom[1],atom[2],(natom,1))
    else:
      atomout=(atom[0],atom[1],atom[2],(natom,2))
#   print atomout[3]
#   sys.stdout.flush()
    latomout[i]=atomout
  return latomout

# BEGIN GAUSSIAN INPUT PROCESSOR SECTION
# The definitions below can extract geometry from most Gaussian inputs
# (including most Z-matrix inputs). 
def is_gaussian_input(li):
  """Returns True if the supplied list of lines looks like a Gaussian input
  and False otherwise. Pretty simple and may be fooled easily."""
  for s in li:
    s=s.strip()
    if (len(s)>0) and s[0]!="%" and s[0]!="#":
      return False
    if s[0]=="#":
      return True
  return False

def find_atom_type(s):
  """Returns the atomic number of the atom specified by the string s, or
  0 for a dummy atom (X) or a ghost atom (Bq), or -100 for an unknown atom."""
  s=s.upper()
  if s.find("-BQ")>=0:
    return 0
  s=(s.split("("))[0]
  s=(s.split("-"))[0]
  if s.isdigit():
    return int(s)
  else:
    while True:
      if s[-1].isdigit():
        s=s[:-1]
      else:
        break
#s has been "standardized", now just look it up:
    if s=="H":
      return 1
    elif s=="HE":
      return 2
    elif s=="LI":
      return 3
    elif s=="BE":
      return 4
    elif s=="B":
      return 5
    elif s=="C":
      return 6
    elif s=="N":
      return 7
    elif s=="O":
      return 8
    elif s=="F":
      return 9
    elif s=="NE":
      return 10
    elif s=="NA":
      return 11
    elif s=="MG":
      return 12
    elif s=="AL":
      return 13
    elif s=="SI":
      return 14
    elif s=="P":
      return 15
    elif s=="S":
      return 16
    elif s=="CL":
      return 17
    elif s=="AR":
      return 18
    elif s=="X" or s=="XX":
      return 0
    else:
      return -100

def resolve_symbol(dlookup,sym):
  """Tries to resolve the expression sym using the lookup dictionary
  lookup."""
#first try to split sym into names,numbers, and arithmetic symbols
  sym=sym.upper()
  unresolved=0
  lsym=[]
  ssym=""
  for schar in sym:
    if schar=="+" or schar=="-" or schar=="*" or schar=="/":
      if len(ssym)>0:
        lsym.append(ssym)
      lsym.append(schar)
      ssym=""
    else:
      ssym+=schar
  lsym.append(ssym)
  for i in range(len(lsym)):
    if dlookup.has_key(lsym[i]):
      lsym[i]=dlookup[lsym[i]] 
  sout="".join(lsym)
  try:
    aout=eval(sout)
  except:
    unresolved=1
    aout=0.0
  sx=str(aout)
  return(sx,unresolved)

def read_gaussian_geometry(li):
  """Extracts the geometry section from the Gaussian input (given as a list
  of lines), removes the (charge,spin) line and other unnecessary data.
  If fragments are specified, includes the fragment information as well."""
  lgeom=[]
  lvar=[]
  lvar2=[]
  ipart=0
  for s in li:
    s=s.strip()
    if ipart==6:
      break
    if ipart==5:
      if len(s)==0:
        ipart=6
      else:
        lvar2.append((s.replace(","," ")).replace("="," "))
    if ipart==4:
      if len(s)==0 or s.upper()=="VARIABLES:" or s.upper()=="CONSTANTS:":
        ipart=5
      else:
        lvar.append((s.replace(","," ")).replace("="," "))
    if ipart==3:
      if len(s)==0 or s.upper()=="VARIABLES:" or s.upper()=="CONSTANTS:":
        ipart=4
      else:
        lgeom.append(s)
    if ipart==2 and len(s)==0:
      ipart=3
    if ipart==1 and len(s)==0:
      ipart=2
    if ipart==0 and len(s)>0 and s[0]=="#":
      ipart=1
  lgeom=lgeom[1:]
# first substitute spaces for commas in lgeom (but not within parentheses)
  ltemp=[]
  for s in lgeom:
    inparentheses=0
    sout=""
    for schar in s:
      if schar=="(":
        inparentheses=1
      if schar==")":
        inparentheses=0
      if schar=="," and inparentheses==0:
        sout+=" "
      else:
        sout+=schar
    ltemp.append(sout)
  lgeom=ltemp
# now substitute line numbers for (possible) atom types in the geometry
# specification, and substitute values for variables using a lookup
# dictionary dlookup:
  dlookup={}
  indeks=1
  for s in lgeom:
    k=s.split()
    s0=k[0]
    s0=(s0.split("("))[0]
    s0=(s0.split("-"))[0]
    s0=s0.upper()
    dlookup[s0]=str(indeks)
    indeks+=1
# check if some variable symbols need resolving
  unresolved=0
  for s in lgeom:
    k=s.split()
    for s1 in k[1:]:
      (sdum,unr)=resolve_symbol(dlookup,s1)
      if unr:
        unresolved=1
  if unresolved:
    for s in lvar:
      k=s.split()
      dlookup[k[0].upper()]=k[1]
# check one last time, sometimes the variable declarations come in two pieces
  unresolved=0
  for s in lgeom:
    k=s.split()
    for s1 in k[1:]:
      (sdum,unr)=resolve_symbol(dlookup,s1)
      if unr:
        unresolved=1
  if unresolved:
    for s in lvar2:
      k=s.split()
      dlookup[k[0].upper()]=k[1]
# lookup dictionary ready
  ltemp=[]
  for s in lgeom:
    k=s.split()
    kout=[k[0]]
    for s1 in k[1:]:
      (sres,unr)=resolve_symbol(dlookup,s1)
      if unr:
        print ("Unresolved symbol",s1,"in geometry specification")
        raise "Error!"
      kout.append(sres)
    sout=" ".join(kout)
    ltemp.append(sout)
  lgeom=ltemp
# symbols resolved, now check for fragment information
  lfrag=[]
  aretherefrag=0
  for s in lgeom:
    k=s.split()
    s0=(k[0]).upper()
    indeks=s0.find("FRAGMENT")
    if indeks>0:
      snum=s0[indeks+9:]
      snum=(snum.split(")"))[0]
      snum=(snum,split(","))[0]
      lfrag.append(int(snum))
      foundfrag=1
    else:
      lfrag.append(0)
      foundfrag=-1
    if foundfrag*aretherefrag==-1:
      print ("Inconsistent fragment specification")
      raise "Error!"
    aretherefrag=foundfrag
# now translate atom types to atomic numbers
  ltemp=[]
  for s in lgeom:
    k=s.split()
    atomtype=find_atom_type(k[0])
    if atomtype==-100:
      print ("Unsupported or unrecognized atom type:",k[0])
      raise "Error!"
    k[0]=str(atomtype)
    sout=" ".join(k)  
    ltemp.append(sout)
  return(ltemp,lfrag)

def gaussian_units_are_au(li):
  """Returns True if the distances in the Gaussian geometry specification are
  in AU and False if they are in Angstrom."""
  sroute=""
  iroute=0
  for s in li:
    s=s.strip()
    if (len(s)>0) and s[0]=="#":
      sroute=s
      iroute=1
      continue
    if iroute==1:
      if len(s)>0:
        sroute=sroute+" "+s
      else:
        break
  sroute=sroute.upper()
  indeks=sroute.find("UNITS")
  if indeks<0:
    return False
  else:
    sroute=sroute[indeks+6:]
    if sroute[0]=="(":
      sroute=(sroute.split(")"))[0]
      sroute=sroute[1:]
    else:
      sroute=(sroute.split(","))[0]
      sroute=(sroute.split())[0]
    if sroute.find("RAD")>=0:
      print ("Angles in radians are not supported")
      raise "Error!"
    if sroute.find("AU")>=0:
      return True
    else:
      return False

def gaussian_input_is_suitable(li):
  """Returns True if the Gaussian input is free from several problems that
  will fool the geometry read."""
  sroute=""
  iroute=0
  for s in li:
    s=s.strip()
    if (len(s)>0) and s[0]=="#":
      sroute=s
      iroute=1
      continue
    if iroute==1:
      if len(s)>0:
        sroute=sroute+" "+s
      else:
        break
  sroute=sroute.upper()
  iproblem=0
  if sroute.find("ONIOM")>=0:
    print ("ONIOM input not supported")
    iproblem+=1 
  indeks=sroute.find("GEOM")
  if indeks>=0:
    sroute=sroute[indeks+5:]
    if sroute[0]=="(":
      sroute=(sroute.split(")"))[0]
      sroute=sroute[1:]
    else:
      sroute=(sroute.split(","))[0]
      sroute=(sroute.split())[0]
    if sroute.find("CHECK")>=0:
      print ("Geom=CheckPoint or Geom=AllCheck not supported")
      iproblem+=1
    if sroute.find("MODELA")>=0:
      print ("Geom=Modela not supported")
      iproblem+=1
  if iproblem==0:
    return True
  else:
    return False 

def vec_norm(vec):
  """Returns the length of a three-dimensional vector."""
  return math.sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])

def scalar_product(vec1,vec2):
  """Returns the scalar product of two 3D vectors."""
  return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]

def vector_product(vec1,vec2):
  """Returns the vector product of two 3D vectors."""
  return [vec1[1]*vec2[2]-vec1[2]*vec2[1],
          vec1[2]*vec2[0]-vec1[0]*vec2[2],
          vec1[0]*vec2[1]-vec1[1]*vec2[0]]

def process_gaussian_geometry(li):
  """Reads in the geometry from the Gaussian input in the 'standardized' form
  (processed by the function read_gaussian_geometry)."""
  lcart=[]
  natoms=0
  for s in li:
    k=s.split()
    if len(k)==4:
#Cartesians
      lcart.append(k)
      natoms+=1
    elif len(k)==1:
#first atom in Z-matrix
      if natoms!=0:
        print ("Incorrect Z-matrix, line ",s)
        raise "Error!"
      else:
        k.extend(["0.0","0.0","0.0"])
        lcart.append(k)
        natoms+=1
    elif k[1]=="0":
#Cartesians mixed into a Z-matrix. If they are really mixed, the results may
#be wrong!
      del k[1]
      lcart.append(k)
      natoms+=1
    elif len(k)==3:
#second atom in Z-matrix
      if natoms!=1 or k[1]!="1":
        print ("Incorrect Z-matrix, line ",s)
        raise "Error!"
      else:
        kout=[k[0],"0.0","0.0",k[2]]
        lcart.append(kout)
        natoms+=1
    elif len(k)==5:
#third atom in Z-matrix
      if natoms!=2 or (k[1]!="1" and k[1]!="2"):
        print ("Incorrect Z-matrix, line ",s)
        raise "Error!"
      else:
        r=float(k[2])
        alpha=float(k[4])
        x=r*math.sin(alpha*math.pi/180.0)
        z=r*math.cos(alpha*math.pi/180.0)
        if k[1]=="1":
          if k[3]!="2":
            print ("Incorrect Z-matrix, line ",s)
            raise "Error!"
          kout=[k[0],str(x),"0.0",str(z)]
          lcart.append(kout)
          natoms+=1 
        else:
          if k[3]!="1":
            print ("Incorrect Z-matrix, line ",s)
            raise "Error!"
          zstart=float(lcart[1][3])
          kout=[k[0],str(x),"0.0",str(zstart-z)]
          lcart.append(kout)
          natoms+=1 
    elif (len(k)==7) or (len(k)==8 and k[7]=="0"):
#general Z-matrix line with R,bond angle,dihedral angle
      if natoms<3 or k[1]==k[3] or k[1]==k[5] or k[3]==k[5]:
        print ("Incorrect Z-matrix, line ",s)
        raise "Error!"
      r=float(k[2])
      alpha=float(k[4])
      phi=float(k[6])
#fetch the Cartesian coordinates of the atoms specified in the line
      nrat1=int(k[1])-1
      nrat2=int(k[3])-1
      nrat3=int(k[5])-1
      cat1=[float(lcart[nrat1][1]),float(lcart[nrat1][2]),
            float(lcart[nrat1][3])]
      cat2=[float(lcart[nrat2][1]),float(lcart[nrat2][2]),
            float(lcart[nrat2][3])]
      cat3=[float(lcart[nrat3][1]),float(lcart[nrat3][2]),
            float(lcart[nrat3][3])]
#find the directions from 1 to 2 and from 2 to 3
      dir12=[(cat2[i]-cat1[i]) for i in range(3)]
      norm12=vec_norm(dir12)
      dir12n=[x/norm12 for x in dir12]
      dir23=[(cat3[i]-cat2[i]) for i in range(3)]
      norm23=vec_norm(dir23)
      dir23n=[x/norm23 for x in dir23]
#orthogonalize dir23n to dir12n
      cos23_12=scalar_product(dir12n,dir23n)
      dir23o12=[(x2-cos23_12*x1) for (x1,x2) in zip(dir12n,dir23n)]
      norm23o12=vec_norm(dir23o12)
      dir23o12n=[x/norm23o12 for x in dir23o12]
#find third direction complementing (dir12n,dir23o12n) to an orthonormal basis
      dircompln=vector_product(dir12n,dir23o12n)
#now we have a nice basis to express the coefficients of the vector
#from atom 1 to the new atom
      x10=r*math.cos(alpha*math.pi/180.0)
      y10=r*math.sin(alpha*math.pi/180.0)*math.cos(phi*math.pi/180.0)
      z10=r*math.sin(alpha*math.pi/180.0)*math.sin(phi*math.pi/180.0)
      dir10=[x10*dir12n[0]+y10*dir23o12n[0]+z10*dircompln[0],
             x10*dir12n[1]+y10*dir23o12n[1]+z10*dircompln[1],
             x10*dir12n[2]+y10*dir23o12n[2]+z10*dircompln[2]]
      cat0=[(dir10[i]+cat1[i]) for i in range(3)]
      kout=[k[0],str(cat0[0]),str(cat0[1]),str(cat0[2])]
      lcart.append(kout)
      natoms+=1 
    elif len(k)==8 and k[7]=="1":
      print ("NYI")
      raise "Error!"
#general Z-matrix line with R,bond angle,second bond angle
      if natoms<3 or k[1]==k[3] or k[1]==k[5] or k[3]==k[5]:
        print ("Incorrect Z-matrix, line ",s)
        raise "Error!"
      r=float(k[2])
      alpha=float(k[4])
      beta=float(k[6])
#fetch the Cartesian coordinates of the atoms specified in the line
      nrat1=int(k[1])-1
      nrat2=int(k[3])-1
      nrat3=int(k[5])-1
      cat1=[float(lcart[nrat1][1]),float(lcart[nrat1][2]),
            float(lcart[nrat1][3])]
      cat2=[float(lcart[nrat2][1]),float(lcart[nrat2][2]),
            float(lcart[nrat2][3])]
      cat3=[float(lcart[nrat3][1]),float(lcart[nrat3][2]),
            float(lcart[nrat3][3])]
    else:
      print ("Incorrect Z-matrix, line ",s)
      raise "Error!"
  return lcart 

def convert_to_bohr(lin):
  """Converts Cartesian coordinates from Angstrom to bohr."""
  convfact=0.529177209
  lout=[]
  for atom in lin:
    atom[1]=str(float(atom[1])/convfact)
    atom[2]=str(float(atom[2])/convfact)
    atom[3]=str(float(atom[3])/convfact)
    lout.append(atom)
  return lout

def omit_dummy_atoms(lcart):
  """Removes dummy atoms from the list of Cartesian coordinates."""
  lpure=[]
  for atom in lcart:
    if atom[0]!="0":
      lpure.append(atom)
  return lpure

def separate_monomers(lin,lfrag,version):
  """Decides which atoms belong to A and which to B."""
  fraga=lfrag[0]
  laa=[]
  lab=[]
  if fraga!=0:
    print ("Using fragments defined in the Gaussian input file.")
    for (atom,fragment) in zip(lin,lfrag):
      attrans=(float(atom[1]),float(atom[2]),float(atom[3]),int(atom[0]))
      if fragment==fraga and atom[0]!="0":
        laa.append(attrans)
      if fragment!=fraga and atom[0]!="0":
        lab.append(attrans)
  else:
    print ("Fragments A and B not specified: attempting to guess.")
    lnondummy=omit_dummy_atoms(lin)
    bestsep=0.0
    bestind=0
    secondbestsep=0.0
    for indsep in range(1,len(lnondummy)):
      mindist=1000.0
      for ata in lnondummy[:indsep]:
        xa=float(ata[1])
        ya=float(ata[2])
        za=float(ata[3])
        for atb in lnondummy[indsep:]:
          xb=float(atb[1])
          yb=float(atb[2])
          zb=float(atb[3])
          rab=math.sqrt((xa-xb)*(xa-xb)+(ya-yb)*(ya-yb)+(za-zb)*(za-zb)) 
          if rab<mindist:
            mindist=rab
      if mindist>bestsep:
        secondbestsep=bestsep
        bestsep=mindist
        bestind=indsep
      elif mindist>secondbestsep:
        secondbestsep=mindist 
    if bestsep-secondbestsep<0.5:
      print ("Could not reliably identify monomers A and B")
      raise "Error!"
    else:
      for atom in lnondummy[:bestind]:
        attrans=(float(atom[1]),float(atom[2]),float(atom[3]),int(atom[0]))
        laa.append(attrans)
      for atom in lnondummy[bestind:]:
        attrans=(float(atom[1]),float(atom[2]),float(atom[3]),int(atom[0]))
        lab.append(attrans)
  ltemp=laa
  laa=analyze_hydrogens(ltemp)
  ltemp=lab
  lab=analyze_hydrogens(ltemp)
  if version==2:
    ltemp=laa
    laa=analyze_atoms2(ltemp,3)
    ltemp=lab
    lab=analyze_atoms2(ltemp,3)
    ltemp=laa
    laa=analyze_atoms(ltemp,4)
    ltemp=lab
    lab=analyze_atoms(ltemp,4)
    ltemp=laa
    laa=analyze_atoms2(ltemp,11)
    ltemp=lab
    lab=analyze_atoms2(ltemp,11)
    ltemp=laa
    laa=analyze_atoms(ltemp,12)
    ltemp=lab
    lab=analyze_atoms(ltemp,12)
  return (laa,lab)
# END GAUSSIAN INPUT PROCESSOR SECTION

def fitted_dispersion(laa,lab,dict,dict8,dicd,aswitch,bswitch,onesystem):
  """Calculates the fitted dispersion energy, either within a molecule
     or between two molecules."""
  dispfit=0.0
  if onesystem:
    for x in laa:
      for y in laa:
        r=math.sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2])) 
        if r>0.0:
          dispfit+=disp_total(r,x[3],y[3],dict,dict8,dicd,aswitch,bswitch)
    dispfit*=0.5
  else:
    for x in laa:
      for y in lab:
        r=math.sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2])) 
        dispfit+=disp_total(r,x[3],y[3],dict,dict8,dicd,aswitch,bswitch)
  return dispfit

def fitted_disp_simple(x,listy,dict,dict8,dicd,aswitch,bswitch):
  """A simplified version of fitted_dispersion where x is just a single atom."""
  dispfit=0.0
  for y in listy:
    r=math.sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]))
    dispfit+=disp_total(r,x[3],y[3],dict,dict8,dicd,aswitch,bswitch)
  return dispfit

def disp_grad_numerical(laa,lab,dict,dict8,dicd,aswitch,bswitch,onesystem,
  atomnumber):
  """Calculates numerically the gradient of the fitted dispersion energy
     with respect to the coordinates of atom (atomnumber)."""
  h=0.001
  h2=h*0.5
  gradient=3*[0.0]
  if onesystem:
    x=list(laa[atomnumber])
    listy=laa[:atomnumber]+laa[atomnumber+1:]
  else:
    if atomnumber>=len(laa):
      x=list(lab[atomnumber-len(laa)])
      listy=laa
    else:
      x=list(laa[atomnumber])
      listy=lab
  for i in range(3):
    x[i]+=h2
    dispa=fitted_disp_simple(x,listy,dict,dict8,dicd,aswitch,bswitch)
    x[i]-=h
    dispb=fitted_disp_simple(x,listy,dict,dict8,dicd,aswitch,bswitch)
    x[i]+=h2
    gradient[i]=(dispa-dispb)/h
  return gradient

def derdamp_TT(r,at1,at2,dicd):
  """Returns the first derivative of the Tang-Toennies damping function f_6
     at separation r."""
# it is assumed below that the damping factor is the geometric mean of atomic
# damping factors:
  alpha=math.sqrt(abs(dicd[at1]*dicd[at2]))
  br=alpha*r
  sum=1.0
  term=1.0
  for i in [1,2,3,4,5,6]:
    term*=br/i
    sum+=term
  dsum=alpha*(sum-term)
  dd=-math.exp(-br)*dsum+alpha*math.exp(-br)*sum
  return dd

def derdamp_TT8(r,at1,at2,dicd):
  """Returns the first derivative of the Tang-Toennies damping function f_8 
     at separation r."""
# it is assumed below that the damping factor is the geometric mean of atomic
# damping factors:
  alpha=math.sqrt(abs(dicd[at1]*dicd[at2]))
  br=alpha*r
  sum=1.0
  term=1.0
  for i in [1,2,3,4,5,6,7,8]:
    term*=br/i
    sum+=term
  dsum=alpha*(sum-term)
  dd=-math.exp(-br)*dsum+alpha*math.exp(-br)*sum
  return dd

def derswitching_function(r,at1,at2,aswitch,bswitch):
  """Returns the first derivative of the dispersion switching (Fermi) function 
  for atoms at1,at2 at separation r."""
  r0=aswitch*(covalent_radius(at1)+covalent_radius(at2))
  dswf=bswitch*math.exp(-bswitch*(r/r0-1.0))/(r0*
    (1.0+math.exp(-bswitch*(r/r0-1.0)))*(1.0+math.exp(-bswitch*(r/r0-1.0))))
  return dswf

def disp_grad_analytical(laa,lab,dict,dict8,dicd,aswitch,bswitch,onesystem,
  atomnumber):
  """Calculates analytically the gradient of the fitted dispersion energy
     with respect to the coordinates of atom (atomnumber)."""
  gradient=3*[0.0]
  if onesystem:
    x=list(laa[atomnumber])
    listy=laa[:atomnumber]+laa[atomnumber+1:]
  else:
    if atomnumber>=len(laa):
      x=list(lab[atomnumber-len(laa)])
      listy=laa
    else:
      x=list(laa[atomnumber])
      listy=lab
  for y in listy:
    r=math.sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]))
    ec6=catomatom(x[3],y[3],dict)/(r*r*r*r*r*r)
    ec8=catomatom8(x[3],y[3],dict8)/(r*r*r*r*r*r*r*r)
    ed6=damp_TT(r,x[3],y[3],dicd)
    ed8=damp_TT8(r,x[3],y[3],dicd)
    esw=switching_function(r,x[3],y[3],aswitch,bswitch)
    dec6=-6.0*catomatom(x[3],y[3],dict)/(r*r*r*r*r*r*r)
    dec8=-8.0*catomatom8(x[3],y[3],dict8)/(r*r*r*r*r*r*r*r*r)
    ded6=derdamp_TT(r,x[3],y[3],dicd)
    ded8=derdamp_TT8(r,x[3],y[3],dicd)
    desw=derswitching_function(r,x[3],y[3],aswitch,bswitch)
    dedr=-esw*(dec6*ed6+ec6*ded6+dec8*ed8+ec8*ded8)-desw*(ec6*ed6+ec8*ed8)
    for i in range(3):
      gradient[i]+=dedr*(x[i]-y[i])/r
  return gradient

#main program begins
#optimized parameters below
nnatoms=[(1,1),(1,3),(1,4),(1,5),(1,6),(1,7),(1,8),(1,9),(1,11),(1,12),(1,13),
         (1,14),(1,15),(1,16),(1,17),2,(3,0),(3,1),(3,2),(4,0),(4,1),(4,2),5,6,
         7,8,9,10,(11,0),(11,1),(11,2),(12,0),(12,1),(12,2),13,14,15,16,17,18]
nnc6=[0.2044,6.1438,0.4117,0.3242,0.0794,0.1792,0.1863,0.0939,5.3056,0.7772,
      0.6309,0.5608,0.1477,0.1272,0.1965,0.0692,39.7169,0.0727,0.5970,12.2461, 
      5.1374,1.2402,0.4142,1.5763,1.0182,0.6889,0.5181,0.2459,45.7198,0.0000,
      4.5841,27.5251,25.8102,2.7405,1.7066,1.1809,3.2334,7.9959,
      5.1041,3.0519]
nnc8=[0.0097,0.0000,0.1598,0.0729,0.0132,0.0040,0.0017,0.0000,0.0013,0.6239,
      0.1697,0.1011,0.5304,0.0674,0.0009,0.0045,37.6958,0.0000,0.0000,4.1540,
      3.8813,0.0408,0.0020,0.1613,0.1481,0.0445,0.1202,0.0312,64.1860,21.2101,
      0.0000,18.3248,34.2080,0.0458,0.1561,0.0000,0.8047,0.3412,
      1.0469,0.6312]
nnbeta=[1.9610,1.6645,1.5344,1.6810,1.4897,1.5787,1.7128,2.2427,1.5209,1.1986,
        1.7140,1.7657,0.5693,1.2785,2.0054,2.3898,0.7598,2.1752,1.2077,1.1788,  
        0.6077,1.7013,2.4917,1.9509,2.3045,2.3927,1.9017,2.1774,0.7055,0.8651,  
        1.2101,0.9876,0.6883,1.7237,1.1894,1.0390,4.1080,2.5881,
        1.5408,1.8462]
#parameters for the JPCL2010 version below
nnatoms0=[(1,6),(1,7),(1,8),(1,9),(1,16),(1,17),2,6,7,8,9,10,16,17,18]
nnc60=[0.0805,0.1851,0.1602,0.0960,0.1740,0.2335,0.0654,1.5290,1.0508,0.8435,
  0.5202,0.2355,8.0173,4.9415,2.9516]
nnc80=[0.0142,0.0026,0.0027,0.0000,0.0572,0.0003,0.0053,0.1649,0.1284,0.0210,
  0.1108,0.0353,0.2105,1.0257,0.6782]
nnbeta0=[1.4632,1.7192,1.7726,2.1695,1.3654,2.0669,2.3164,1.9509,2.2992,
  2.4256,1.9457,2.1363,2.7729,1.5878,1.8374]
#optimized parameters above
dict={}
dict8={}
dicd={}
if (len(sys.argv)<3):
  print ("Usage: python disp4dldf.py dir_with_saptdft_outputs version                \n   or: python disp4dldf.py input_file version optional_arguments")
  print ("Version: 1 for Podeszwa et al., JPCL 1, 550 (2010)")
  print ("         2 for Patkowski et al., to be published")
  print ("If the optional third argument is set to 1, the intramolecular")
  print ("dispersion energy for a single system is computed. Normally,")
  print ("the intermolecular dispersion energy is computed for a dimer.")
  print ("The optional arguments 'a' and 'n' request analytical and")
  print ("numerical gradients, respectively.")
  raise "Error!"
if (not(os.path.exists(sys.argv[1]))):
  print ("Input file/directory "+sys.argv[1]+" not found!")
  raise "Error!"
onesystem=False
angrad=False
numgrad=False
if (os.path.isdir(sys.argv[1])):
  mode=1
  if len(sys.argv)>3:
    print ("The optional arguments are not compatible with a directory read.")
    raise "Error!"
else:
  mode=2
  if len(sys.argv)>3:
    for sarg in sys.argv[3:]:
      if sarg=="1":
        onesystem=True
      elif sarg=="a":
        angrad=True
      elif sarg=="n":
        numgrad=True
      else:
        print ("Optional argument",sarg,"not recognized.")
        raise "Error!"
version=int(sys.argv[2])
if version!=1 and version!=2:
  print ("Select dispersion version 1 or 2.")
  raise "Error!"
if version==1:
  for i in range(len(nnatoms0)):
    dict[nnatoms0[i]]=nnc60[i]
    dict8[nnatoms0[i]]=nnc80[i]
    dicd[nnatoms0[i]]=nnbeta0[i]
else:
  for i in range(len(nnatoms)):
    dict[nnatoms[i]]=nnc6[i]
    dict8[nnatoms[i]]=nnc8[i]
    dicd[nnatoms[i]]=nnbeta[i]
#default values for the switching function parameters:
aswitch=1.0
bswitch=8.8
if (os.path.exists("switch.par")):
  lsw=[line.rstrip() for line in open("switch.par","r")]
  aswitch=float(lsw[0])
  bswitch=float(lsw[1])
rmse=0.0
mues=0.0
mures=0.0
nsyst=0
icountf=0
latoma=[]
latomb=[]
results=[]
if (mode==1):
  f=os.popen("ls -1 "+sys.argv[1]+"/*.out")
  li=f.readlines()
  f.close()
  for s in li:
    s=s.strip("\n")
    f1=open(s,"r")
    li1=f1.readlines()
    f1.close()
#   print "Reading file",s
    (laa,lab,res)=read_saptdft_output(li1,version)
    latoma.append(laa)
    latomb.append(lab)
    results.append(res)
    icountf=icountf+1
else:
  f=open(sys.argv[1],"r")
  li=f.readlines()
  f.close()
  if is_gaussian_input(li):
    print( sys.argv[1],"appears to be a Gaussian input file.")
    if not(gaussian_input_is_suitable(li)):
      print ("Input file could not be processed.")
      raise "Error!"
    (lgeom,lfrag)=read_gaussian_geometry(li)
    lcart=process_gaussian_geometry(lgeom)
    if not(gaussian_units_are_au(li)):
      lin=convert_to_bohr(lcart)
    else:
      lin=lcart
    if not(onesystem):
      (laa,lab)=separate_monomers(lin,lfrag,version)
    else:
      laa=[]
      for atom in lin:
        attrans=(float(atom[1]),float(atom[2]),float(atom[3]),int(atom[0]))
        if atom[0]!="0":
          laa.append(attrans)
      ltemp=laa
      laa=analyze_hydrogens(ltemp)
      if version==2:
        ltemp=laa
        laa=analyze_atoms2(ltemp,3)
        ltemp=laa
        laa=analyze_atoms(ltemp,4)
        ltemp=laa
        laa=analyze_atoms2(ltemp,11)
        ltemp=laa
        laa=analyze_atoms(ltemp,12)
      lab=[]
  else:
#    print (sys.argv[1],"is not a Gaussian input. Assuming the dimer.cnf format.")
    (laa,lab)=read_geometry(li,version)
    if onesystem:
      laa.extend(lab)
  latoma.append(laa)
  latomb.append(lab)
  icountf=icountf+1
for i in range(icountf):
  dispfit=fitted_dispersion(latoma[i],latomb[i],dict,dict8,dicd,
    aswitch,bswitch,onesystem)
  if (numgrad):
    nrat=len(latoma[i])+len(latomb[i])
    print ("Numerical gradient:")
    for atomnumber in range(nrat):
      print (atomnumber,disp_grad_numerical(latoma[i],latomb[i],dict,dict8,
        dicd,aswitch,bswitch,onesystem,atomnumber))
    print ("End numerical gradient")
  if (angrad):
    nrat=len(latoma[i])+len(latomb[i])
    print ("Analytical gradient:")
    for atomnumber in range(nrat):
      print (atomnumber,disp_grad_analytical(latoma[i],latomb[i],dict,dict8,
        dicd,aswitch,bswitch,onesystem,atomnumber))
    print ("End analytical gradient")
  if (mode==1):
    print ("==============================================================")
    print ("Working with file",li[i])
    print ("Benchmark E(disp+ex-disp) in kcal/mol: ",results[i])
    rmsec=(dispfit-results[i])*(dispfit-results[i])/(results[i]*results[i])
    rmse=rmse+rmsec
    errorx=100.0*(dispfit-results[i])/(results[i])
    mures=mures+abs(errorx)
    mues=mues+abs(dispfit-results[i])
    serrorx="%6.2f" % errorx
  nsyst=nsyst+1
  if version==1:
    print ("Fitted dispersion energy (JPCL2010):       ",dispfit)
  else:
    print (round(dispfit,8))
  if (mode==1):
    print ("Percent error:                             ",serrorx)
if (mode==1):
  mue=mues/nsyst
  mure=mures/nsyst
  smue="%6.2f" % mue
  smure="%6.2f" % mure
  print ("==============================================================")
  print ("==============================================================")
  print ("Total RMSE:                            ",rmse)
  print ("Dispersion MUE (kcal/mol):             ",smue)
  print ("Dispersion MURE (%):                   ",smure)
