%Chk=ArHFB_1
%Mem=6GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
 H1        0.00000000000          0.00000000000          2.30429326844
 F1        0.00000000000          0.00000000000          3.22125160551
 Ar-Bq         0.00000000000          0.00000000000          0.00000000000

@aug-cc-pVTZ.gbs/N

