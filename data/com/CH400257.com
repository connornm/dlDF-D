%Chk=CH400257
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.64422640        0.00000000       -0.90053223
  1        -0.99534192        0.41059314       -0.25825123
  1        -0.11364261       -1.03649085        0.37250240
  1         0.46475812        0.62589771        0.78628106
  0        -0.42629406        0.00000000        0.59589539
  0         0.65863235       -0.27169550        0.17088863
  0         0.07519898        0.68586120       -0.24649030
  0        -0.30753727       -0.41416569       -0.52029371
  6         0.00000000        0.00000000        3.25190421
  1         0.16063186        1.09536343        3.23286945
  1        -0.59958960       -0.30046541        2.37088276
  1         0.97875343       -0.51706508        3.22599312
  1        -0.53979568       -0.27783294        4.17787151
  0        -0.10629246       -0.72481805        3.26449979
  0         0.39675724        0.19882237        3.83488903
  0        -0.64765550        0.34214955        3.26904996
  0         0.35719072        0.18384613        2.63917807

@aug-cc-pVTZ.gbs/N

