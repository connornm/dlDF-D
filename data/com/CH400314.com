%Chk=CH400314
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.98671828        0.00000000        0.50236739
  1        -0.51760870        0.95748853        0.20318194
  1         0.14166610       -0.11855577       -1.09172396
  1        -0.61077568       -0.83893276        0.38617463
  0        -0.65292596        0.00000000       -0.33242387
  0         0.34250927       -0.63358421       -0.13444847
  0        -0.09374254        0.07845009        0.72240976
  0         0.40415923        0.55513411       -0.25553742
  6         0.00000000        0.00000000        3.33838631
  1         0.01761834        0.00022787        4.44548861
  1        -1.03221078       -0.19014848        2.98572439
  1         0.34270651        0.98390815        2.96360242
  1         0.67188592       -0.79398754        2.95872982
  0        -0.01165831       -0.00015078        2.60580048
  0         0.68302901        0.12582404        3.57174788
  0        -0.22677393       -0.65106645        3.58638631
  0        -0.44459677        0.52539320        3.58961058

@aug-cc-pVTZ.gbs/N

