#!/usr/bin/env python
# Fix old mistake I made
#Usage:
# 	./rename.py NAME FILETYPE

import os, re, sys
key=sys.argv[1]
ext=sys.argv[2]
fillnum=5

for filename in os.listdir('.'):
	filenew = re.sub(key, '', filename)
	filenew = re.sub('.'+ext, '', filenew)
	filenew = filenew.zfill(fillnum)
	filenew=key+filenew+'.'+ext
	os.system('mv '+filename+' '+filenew)

os.system('mv '+key+'rename.py.'+ext+' rename.py')
