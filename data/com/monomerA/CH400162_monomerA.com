%Chk=CH400162
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.92566775        0.00000000       -0.60755672
  1        -0.03526024       -0.91743856        0.61890965
  1        -0.00923867        0.89005886        0.65855586
  1        -0.88116884        0.02737970       -0.66990879
  6-Bq        0.00000000       0.00000000       3.52017323
  1-Bq        0.81189980      -0.45694685       4.11850779
  1-Bq       -0.95099097      -0.05650693       4.08445073
  1-Bq        0.24477351       1.06010634       3.31463387
  1-Bq       -0.10568234      -0.54665257       2.56310052

@aug-cc-pVTZ.gbs/N

