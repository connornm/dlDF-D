%Chk=BLIND00982
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.35191959        0.92366487        0.25318910
  6         1.10128679        0.28847253       -0.77676533
 16        -0.49877921        0.09144741        1.54634119
 16         1.27274756       -1.46597084       -0.94144080
  6        -0.98215822       -1.39511647        0.70932588
  6        -0.27678114       -2.01101740       -0.27594225
  6         0.45395940        2.29410315        0.12361285
  7         1.72631752        1.06087356       -1.63151282
 16         1.47687933        2.65648020       -1.22810615
  6        -2.20136034       -1.94385223        1.20280830
  7        -3.18597797       -2.36546791        1.65442914
  6        -0.73255253       -3.23327033       -0.85240749
  7        -1.05051793       -4.24415780       -1.32982831
  6        -0.12200255        3.29944073        0.93769678
  7        -0.57922056        4.13774671        1.60048684
  6         0.13789185       -0.57515860        9.85953434
  6         1.25428921        0.24872616        9.54220565
 16        -1.54062186       -0.18005285        9.52041811
 16         1.13262948        1.80961724        8.71552173
  6        -1.35667384        0.74522574        8.01903717
  6        -0.29593370        1.53441598        7.70335830
  6         0.56758806       -1.69797102       10.53773538
  7         2.43808313       -0.16845758        9.91963075
 16         2.28737880       -1.61091466       10.73665135
  6        -2.47864623        0.62734635        7.14819156
  7        -3.42566790        0.51925290        6.48286549
  6        -0.26817408        2.27619131        6.48559778
  7        -0.21592058        2.92027294        5.51949501
  6        -0.20780037       -2.76767136       11.04732580
  7        -0.82948456       -3.64929551       11.48015781

@aug-cc-pVTZ.gbs/N
