%Chk=NH3_2
%Mem=6GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=20)

M05 opt

0 1
 N2        0.00000000000          0.06336956769          3.18174908834
 H4       -0.80956498945         -0.11963738786          3.76165221007
 H5        0.00000000000         -0.64143473095          2.45250209257
 H6        0.80956498945         -0.11963738786          3.76165221007
 N1        0.00000000000         -0.06336956769          0.02547991166
 H1       -0.80956498945          0.11963738786         -0.55442321007
 H2        0.80956498945          0.11963738786         -0.55442321007
 H3        0.00000000000          0.64143473095          0.75472690743

@aug-cc-pVTZ.gbs/N
