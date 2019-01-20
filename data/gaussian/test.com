%Chk=test
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1

@blind-aug.gbs/N

