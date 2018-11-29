%Chk=ArHFB_3
%Mem=6GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
 H1        0.00000000000          0.00000000000          3.36264726844
 F1        0.00000000000          0.00000000000          4.27960560551
 Ar-Bq         0.00000000000          0.00000000000          0.00000000000

@aug-cc-pVTZ.gbs/N

