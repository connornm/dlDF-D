%Chk=BLIND00916
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.37871040        0.73884756       -0.59312024
  6         0.97300950       -0.53504984       -0.81632565
 16        -1.27772406        1.00540545       -0.06996786
 16         0.13821856       -2.07664584       -0.56895106
  6        -1.52711877       -0.42330581        0.95018335
  6        -0.96520467       -1.64452734        0.74890100
  6         1.27817354        1.73402411       -0.91853217
  7         2.20451963       -0.55561271       -1.26445509
 16         2.74369651        1.00343210       -1.48713073
  6        -2.43906973       -0.19381686        2.02105679
  7        -3.20612983        0.03474716        2.86401011
  6        -1.27003267       -2.74406003        1.60450002
  7        -1.51182778       -3.67436212        2.25795183
  6         1.09687471        3.13709912       -0.85787177
  7         0.96245809        4.29114799       -0.82324103
  6        -0.27233240        0.73104427        5.70564025
  6        -0.83855562       -0.54160322        5.99819156
 16         1.28415431        0.99297505        4.93281555
 16        -0.07058826       -2.08548598        5.59664393
  6         1.35381668       -0.41884093        3.86228070
  6         0.81672397       -1.63869796        4.12897205
  6        -1.09809173        1.72774014        6.18524645
  7        -1.98443713       -0.55984943        6.63418108
 16        -2.46417254        0.99940139        6.96510056
  6         2.08858642       -0.17753785        2.66521098
  7         2.71610413        0.05989474        1.71597954
  6         0.97073457       -2.72489985        3.21771571
  7         1.09620540       -3.64507908        2.51882131
  6        -0.91268089        3.13008905        6.12037129
  7        -0.77229826        4.28338625        6.08439603

@aug-cc-pVTZ.gbs/N
