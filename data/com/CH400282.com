%Chk=CH400282
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.09733974        0.00000000        0.14775467
  1        -0.47216047        0.67397902        0.74081221
  1        -0.23470075        0.35340036       -1.02274615
  1        -0.39047852       -1.02737939        0.13417928
  0        -0.72612580        0.00000000       -0.09777143
  0         0.31243551       -0.44598180       -0.49020630
  0         0.15530493       -0.23385020        0.67676613
  0         0.25838537        0.67983201       -0.08878840
  6         0.00000000        0.00000000        4.10070396
  1        -0.07389459        0.61132543        3.18048292
  1        -0.31156013       -1.03921063        3.87944050
  1         1.04649259        0.00248626        4.46241244
  1        -0.66103787        0.42539894        4.88047997
  0         0.04889713       -0.40452300        4.70962772
  0         0.20616391        0.68766091        4.24711723
  0        -0.69247950       -0.00164519        3.86135615
  0         0.43741846       -0.28149272        3.58471473

@aug-cc-pVTZ.gbs/N

