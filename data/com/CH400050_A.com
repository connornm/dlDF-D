%Chk=CH400050
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10718714        0.00000000       -0.01107251
  1        -0.36232839       -0.79768799        0.67705093
  1        -0.36552114        0.98202010        0.35779431
  1        -0.37933761       -0.18433212       -1.02377273
  6-Bq        0.00000000       0.00000000       4.51907569
  1-Bq        0.17134225      -0.98644701       4.04626485
  1-Bq        0.81555396       0.69254720       4.23405627
  1-Bq       -0.96946252       0.41035081       4.17594220
  1-Bq       -0.01743370      -0.11645100       5.62003945

@aug-cc-pVTZ.gbs/N

