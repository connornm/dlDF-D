%Chk=BLIND01144
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.57601287       -0.38943589        0.74676486
  6        -0.53224948        1.03319358        0.74069892
 16         0.72105621       -1.44438671        0.20527876
 16         0.84105754        1.98130878        0.14934983
  6         1.42005407       -0.45631344       -1.09057435
  6         1.46730304        0.90201812       -1.10920153
  6        -1.76284025       -0.81485994        1.30868869
  7        -1.56475410        1.67025298        1.23663193
 16        -2.68456093        0.57119675        1.79253206
  6         1.99470622       -1.23293583       -2.13829675
  7         2.47052279       -1.90842092       -2.95592228
  6         2.09679981        1.60108100       -2.18110957
  7         2.62436663        2.21079034       -3.01827973
  6        -2.20948877       -2.14185977        1.52021875
  7        -2.58797937       -3.22438958        1.70987633
  6         0.08957312       -0.98846385        6.79063873
  6        -1.14368534       -0.41062432        7.20416227
 16         1.01900237       -0.49171046        5.38425544
 16        -1.93389839        0.94001280        6.37595273
  6         0.64869800        1.24201025        5.34948088
  6        -0.52377454        1.80656548        5.74220162
  6         0.40844546       -2.04468138        7.62011336
  7        -1.73908480       -0.94082283        8.24447427
 16        -0.84803954       -2.23272665        8.79939745
  6         1.70768579        2.02959823        4.81157784
  7         2.58905041        2.62647677        4.34453882
  6        -0.73733649        3.21206559        5.62861171
  7        -0.96569351        4.34717547        5.52709936
  6         1.53579522       -2.89953575        7.56116432
  7         2.44998058       -3.61636376        7.52219456

@aug-cc-pVTZ.gbs/N
