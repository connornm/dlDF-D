%Chk=C2H2ClFB_1
%Mem=6GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
 Cl        0.00000000000          0.00000000000          2.31244547927
  F        0.00000000000          0.00000000000          3.96921636506
 C1-Bq         0.00000000000          0.60529304614          0.00100613208
 C2-Bq         0.00000000000         -0.60529304614          0.00100613208
 H1-Bq         0.00000000000          1.67189106733         -0.01198993198
 H2-Bq         0.00000000000         -1.67189106733         -0.01198993198

@aug-cc-pVTZ.gbs/N

