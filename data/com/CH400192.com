%Chk=CH400192
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.09209685        0.00000000        0.18251144
  1        -0.22169433        0.58659629       -0.91254719
  1        -0.35146428       -1.04113005       -0.13604053
  1        -0.51893824        0.45453376        0.86607628
  0        -0.72265651        0.00000000       -0.12077050
  0         0.14669839       -0.38815937        0.60384586
  0         0.23256907        0.68893103        0.09002002
  0         0.34338904       -0.30077166       -0.57309538
  6         0.00000000        0.00000000        4.18844971
  1         0.32308652        0.03612883        3.13000946
  1        -1.00995246       -0.44953276        4.25091792
  1         0.71381564       -0.61360638        4.77149379
  1        -0.02694970        1.02701031        4.60137768
  0        -0.21379109       -0.02390697        4.88883515
  0         0.66830036        0.29746242        4.14711358
  0        -0.47234228        0.40603234        3.80264089
  0         0.01783301       -0.67958778        3.91520922

@aug-cc-pVTZ.gbs/N

