%Chk=CH400101
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.80794463       0.00000000       0.75710727
  1-Bq       -0.93944247      -0.35957083       0.46275546
  1-Bq       -0.14717760       1.02852306      -0.38270750
  1-Bq        0.27867544      -0.66895223      -0.83715523
  6         0.00000000        0.00000000        2.42052830
  1        -0.58598691       -0.32565532        3.30175127
  1        -0.68926721        0.34959229        1.62763367
  1         0.67662361        0.82651088        2.71215513
  1         0.59863051       -0.85044785        2.04057314

@aug-cc-pVTZ.gbs/N

