%Chk=N2_HFB_3
%Mem=6GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
 H1        0.00000000000          0.00000000000          3.34779372478
 F1        0.00000000000          0.00000000000          4.28039454133
 N1-Bq         0.00000000000          0.00000000000         -0.55073063835
 N2-Bq         0.00000000000          0.00000000000          0.55073063835

@aug-cc-pVTZ.gbs/N

