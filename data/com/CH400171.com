%Chk=CH400171
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.08858818        0.00000000        0.20239007
  1        -0.43655206        0.96293392        0.32888680
  1        -0.17358697       -0.13233565       -1.08551407
  1        -0.47844915       -0.83059827        0.55423721
  0        -0.72033476        0.00000000       -0.13392448
  0         0.28887290       -0.63718750       -0.21762922
  0         0.11486504        0.08756844        0.71830058
  0         0.31659682        0.54961906       -0.36674689
  6         0.00000000        0.00000000        2.75542911
  1        -0.63121101        0.36717374        1.92311783
  1         0.90234860        0.63573742        2.84255480
  1        -0.57444794        0.04294956        3.70102461
  1         0.30331035       -1.04586072        2.55501918
  0         0.41768159       -0.24296425        3.30618169
  0        -0.59709731       -0.42067678        2.69777676
  0         0.38012063       -0.02842035        2.12971469
  0        -0.20070491        0.69206138        2.88804329

@aug-cc-pVTZ.gbs/N

