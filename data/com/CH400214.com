%Chk=CH400214
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.59940581        0.00000000       -0.93096650
  1        -0.80125171        0.76030551       -0.07692322
  1        -0.45269489       -0.99964902        0.14749622
  1         0.65454079        0.23934350        0.86039351
  0        -0.39663562        0.00000000        0.61603419
  0         0.53020001       -0.50310531        0.05090123
  0         0.29955485        0.66148242       -0.09760041
  0        -0.43311924       -0.15837711       -0.56933501
  6         0.00000000        0.00000000        4.63153988
  1         0.41174577        0.25699982        5.62672955
  1         0.67505779       -0.71819000        4.12707376
  1        -1.00102573       -0.45731987        4.75316217
  1        -0.08577783        0.91851004        4.01919405
  0        -0.27245821       -0.17006055        3.97300828
  0        -0.44669564        0.47523685        4.96535251
  0         0.66239341        0.30261527        4.55106063
  0         0.05676045       -0.60779157        5.03673810

@aug-cc-pVTZ.gbs/N

