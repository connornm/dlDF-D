%Chk=BLIND02507
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.39256007       -0.21813557       -0.91619996
  6        -0.31040862       -1.30365160       -0.32179554
 16         0.02349405        1.48248181       -0.67082792
 16        -1.66683088       -1.11228763        0.79984916
  6        -0.56072534        1.44699536        1.00295536
  6        -1.23191777        0.41701277        1.58285245
  6         1.36179338       -0.70481240       -1.77003582
  7         0.05633791       -2.51606410       -0.65856391
 16         1.29029005       -2.43677805       -1.77301118
  6        -0.29929886        2.65501065        1.71249389
  7        -0.08482789        3.66745745        2.24178274
  6        -1.70115784        0.50908707        2.92645880
  7        -2.12679463        0.55149882        4.00718110
  6         2.27571909        0.02400898       -2.56924019
  7         3.02682824        0.60606758       -3.23867185
  6        -1.00551715        0.13399598       10.82712805
  6        -0.22787745        1.23199591       10.36290902
 16        -0.49170951       -1.54660009       10.81639649
 16         1.41553963        1.08325052        9.72131455
  6         1.25047652       -1.35243177       11.08324590
  6         2.00425836       -0.30900032       10.64654488
  6        -2.25095479        0.58070302       11.22037242
  7        -0.78537519        2.41797166       10.38566475
 16        -2.34450980        2.29112717       10.95481026
  6         1.83442702       -2.44697129       11.78472627
  7         2.26878643       -3.37239515       12.33800341
  6         3.41160671       -0.27257732       10.87454910
  7         4.56431510       -0.21521564       11.01202322
  6        -3.33438227       -0.17391634       11.73213696
  7        -4.23705952       -0.77883385       12.14486509

@aug-cc-pVTZ.gbs/N

