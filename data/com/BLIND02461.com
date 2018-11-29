%Chk=BLIND02461
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1

@aug-cc-pVTZ.gbs/N

