%Chk=CH400090
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.93248505        0.00000000       -0.59704069
  1         0.06177548        0.78247461        0.78096302
  1        -0.86252482        0.20719739       -0.66265084
  1        -0.13173571       -0.98967200        0.47872851
  6-Bq        0.00000000       0.00000000       3.56609534
  1-Bq       -0.81354013       0.59506128       4.02439708
  1-Bq       -0.43511593      -0.85033755       3.00610772
  1-Bq        0.58014911       0.63993350       2.87334857
  1-Bq        0.66850695      -0.38465723       4.36052799

@aug-cc-pVTZ.gbs/N

