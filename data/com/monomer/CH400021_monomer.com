%Chk=CH400021
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.05761844        0.00000000        0.32776392
  1        -0.37678286       -1.04070071       -0.03102686
  1        -0.60721107        0.59127634        0.71251173
  1        -0.07362452        0.44942437       -1.00924879
  6-Bq        0.00000000       0.00000000       2.55053616
  1-Bq       -0.85337743      -0.16170898       1.86381705
  1-Bq        0.70113076      -0.85366912       2.47536537
  1-Bq        0.52513393       0.93445060       2.27300709
  1-Bq       -0.37288727       0.08092750       3.58995512

@aug-cc-pVTZ.gbs/N

