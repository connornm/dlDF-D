%Chk=CH400192
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.09209685        0.00000000        0.18251144
  1        -0.22169433        0.58659629       -0.91254719
  1        -0.35146428       -1.04113005       -0.13604053
  1        -0.51893824        0.45453376        0.86607628
  6         0.00000000        0.00000000        4.18844971
  1         0.32308652        0.03612883        3.13000946
  1        -1.00995246       -0.44953276        4.25091792
  1         0.71381564       -0.61360638        4.77149379
  1        -0.02694970        1.02701031        4.60137768

@aug-cc-pVTZ.gbs/N

