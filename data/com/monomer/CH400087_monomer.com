%Chk=CH400087
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.90560308        0.00000000        0.63707851
  1        -0.80881860       -0.55987681        0.50826802
  1        -0.32737240        1.04297670       -0.17610472
  1         0.23058792       -0.48309989       -0.96924180
  6-Bq        0.00000000       0.00000000       2.62220064
  1-Bq       -0.28870513      -0.54232178       1.70104789
  1-Bq       -0.36055749       1.04524386       2.56346800
  1-Bq       -0.45335475      -0.49797168       3.50110885
  1-Bq        1.10261737      -0.00495040       2.72317782

@aug-cc-pVTZ.gbs/N

