%Chk=CH400208
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.87071193       0.00000000       0.68399320
  1-Bq       -0.46994338      -1.00256600       0.00076510
  1-Bq       -0.73673965       0.75321559       0.34039230
  1-Bq        0.33597110       0.24935042      -1.02515060
  6         0.00000000        0.00000000        4.07233275
  1         0.44730993        0.87514549        3.56241104
  1        -0.95907617        0.29573313        4.53998730
  1        -0.18217072       -0.80433638        3.33352102
  1         0.69393696       -0.36654224        4.85341164

@aug-cc-pVTZ.gbs/N
