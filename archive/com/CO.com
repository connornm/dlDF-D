%Chk=ArHF_1
%Mem=6GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 2
 H1        0.00000000000          0.00000000000          0.00000000000
 C1        0.00000000000          0.00000000000          1.00000000000

@aug-cc-pVTZ.gbs/N
