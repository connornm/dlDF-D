%Chk=CH400202
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.99665825        0.00000000       -0.48234664
  1        -0.73454870       -0.48662322       -0.67053859
  1         0.05253140       -0.55651484        0.95578117
  1        -0.31464095        1.04313806        0.19710407
  6-Bq        0.00000000       0.00000000       2.86364946
  1-Bq       -0.34558856      -0.24485207       3.88668512
  1-Bq        0.23665345      -0.93527349       2.32028530
  1-Bq       -0.79796373       0.54757310       2.32568706
  1-Bq        0.90689884       0.63255246       2.92194034

@aug-cc-pVTZ.gbs/N

