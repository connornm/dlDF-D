%Chk=CH400179
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.75579035       0.00000000      -0.80917668
  1-Bq        0.20484912      -0.83611842       0.69636828
  1-Bq        0.04885491       0.95935715       0.55066596
  1-Bq       -1.00949439      -0.12323874      -0.43785757
  6         0.00000000        0.00000000        3.34613958
  1        -0.04987840       -0.24495120        2.26748440
  1        -0.49814907       -0.79854070        3.92937728
  1         1.05945646        0.07872117        3.65813470
  1        -0.51142899        0.96477074        3.52956196

@aug-cc-pVTZ.gbs/N

