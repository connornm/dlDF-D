%Chk=CH400010
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.02678609       0.00000000       0.41436250
  1-Bq       -0.07708056       0.76657767      -0.79523784
  1-Bq       -0.22641106      -0.99696055      -0.42519835
  1-Bq       -0.72329447       0.23038288       0.80607369
  6         0.00000000        0.00000000        4.17586654
  1         0.80272155        0.68856370        4.50374474
  1        -0.95567426        0.29418011        4.65139860
  1        -0.10418287        0.05227014        3.07477630
  1         0.25713558       -1.03501395        4.47354654

@aug-cc-pVTZ.gbs/N

