%Chk=CH400064
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.57395925        0.00000000       -0.94686680
  1         0.33143109       -0.84621959        0.63249648
  1         0.17400505        0.95250405        0.53707004
  1        -1.07939539       -0.10628446       -0.22269973
  0        -0.37979725        0.00000000        0.62655565
  0        -0.21931281        0.55995591       -0.41853220
  0        -0.11514170       -0.63028590       -0.35538713
  0         0.71425176        0.07032999        0.14736368
  6         0.00000000        0.00000000        4.73782371
  1         1.01648958       -0.43192864        4.81638926
  1         0.02627468        0.88767993        4.07652201
  1        -0.69024284       -0.75525374        4.31457133
  1        -0.35252142        0.29950245        5.74381222
  0        -0.67262606        0.28581352        4.68583573
  0        -0.01738634       -0.58739083        5.17541674
  0         0.45674381        0.49976248        5.01789601
  0         0.23326859       -0.19818516        4.07214635

@aug-cc-pVTZ.gbs/N

