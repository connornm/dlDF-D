%Chk=BLIND00585
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.43795705       -0.92099982       -0.03258612
  6        -0.33758606       -0.38135774        1.28063082
 16        -0.62233331        0.02290439       -1.50349294
 16        -0.35040611        1.35111067        1.64529121
  6         0.29753680        1.47585775       -1.07118798
  6         0.40234887        2.00082201        0.17824030
  6        -0.44043632       -2.29985020        0.03120041
  7        -0.26491754       -1.22965814        2.27710749
 16        -0.34633824       -2.78606398        1.69239524
  6         0.91010877        2.10481130       -2.19380219
  7         1.37028974        2.59445439       -3.14231479
  6         1.13072886        3.20463755        0.41111938
  7         1.68862977        4.19771678        0.64234730
  6        -0.53541054       -3.22899507       -1.03314829
  7        -0.62177567       -4.00571472       -1.89358048
  6         0.75170607        0.63631047        5.99350786
  6        -0.19327055        1.21084602        5.09754204
 16         0.39168854       -0.60261587        7.18680814
 16        -1.89843967        0.74106384        5.01822438
  6        -0.88327211       -1.50119241        6.34350653
  6        -1.79078531       -0.96522051        5.48513716
  6         1.97730794        1.25038600        5.83190804
  7         0.22856688        2.16958274        4.30967333
 16         1.83389335        2.48167706        4.62022279
  6        -0.91239938       -2.88448990        6.68516325
  7        -0.91156582       -3.99973796        7.01290412
  6        -2.80884405       -1.77057901        4.89446296
  7        -3.67219980       -2.38217109        4.41332047
  6         3.18583826        0.99686605        6.52498132
  7         4.18394137        0.80856031        7.09016008

@aug-cc-pVTZ.gbs/N

