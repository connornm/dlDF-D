%Chk=CH400029
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.01125962        0.00000000        0.45093230
  1        -0.48237591       -0.98106778        0.17551480
  1        -0.61045950        0.79948935        0.46275473
  1         0.08157579        0.18157843       -1.08920183
  0        -0.66916532        0.00000000       -0.29838852
  0         0.31919521        0.64918694       -0.11614072
  0         0.40395000       -0.52903383       -0.30621159
  0        -0.05397990       -0.12015311        0.72074083
  6         0.00000000        0.00000000        4.24641838
  1        -0.99287452       -0.01698929        3.75662505
  1         0.53524225        0.92710386        3.96361083
  1        -0.13010601       -0.03039067        5.34557023
  1         0.58773828       -0.87972390        3.91986740
  0         0.65699964        0.01124206        4.57052181
  0        -0.35417765       -0.61347823        4.43355628
  0         0.08609305        0.02010995        3.51909348
  0        -0.38891504        0.58212621        4.46250195

@aug-cc-pVTZ.gbs/N

