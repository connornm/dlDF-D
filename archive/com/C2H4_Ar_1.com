%Chk=C2H4_Ar_1
%Mem=6GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
 Ar        0.00000000000          0.00000000000          3.30000000000
 C1        0.66807781003          0.00000000000          0.00000000000
 C2       -0.66807781003          0.00000000000          0.00000000000
 H1        1.23491924972          0.92453298021          0.00000000000
 H2        1.23491924972         -0.92453298021          0.00000000000
 H3       -1.23491924972          0.92453298021          0.00000000000
 H4       -1.23491924972         -0.92453298021          0.00000000000

@aug-cc-pVTZ.gbs/N

