%Chk=ArHFB_2
%Mem=6GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
 H1        0.00000000000          0.00000000000          2.56888226844
 F1        0.00000000000          0.00000000000          3.48584060551
 Ar-Bq         0.00000000000          0.00000000000          0.00000000000

@aug-cc-pVTZ.gbs/N

