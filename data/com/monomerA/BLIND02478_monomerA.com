%Chk=BLIND02478
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=130)

M05 opt

0 1
  6         0.87127980        0.08794027       -0.52368590
  6         0.59648787       -1.23923006       -0.08901940
 16         0.37859721        1.55072984        0.31656834
 16        -0.32519304       -1.63486304        1.36994159
  6        -1.15405032        1.01288697        1.02779171
  6        -1.42961061       -0.25084137        1.44576451
  6         1.64444373        0.05127414       -1.66656663
  7         1.08948654       -2.22167610       -0.80289063
 16         1.97326053       -1.60214535       -2.07018779
  6        -2.10261243        2.06627365        1.17508759
  7        -2.83337713        2.96203933        1.29665989
  6        -2.68182412       -0.56874561        2.04978451
  7        -3.67605421       -0.86626560        2.57315294
  6         2.15158208        1.13761184       -2.42012455
  7         2.58370418        2.01713961       -3.04526238
  6-Bq       -0.87127980      -0.08794023       9.63827511
  6-Bq       -0.59648780       1.23923009       9.20360862
 16-Bq       -0.37859730      -1.55072981       8.79802085
 16-Bq        0.32519310       1.63486304       7.74464762
  6-Bq        1.15405024      -1.01288702       8.08679745
  6-Bq        1.42961060       0.25084131       7.66882466
  6-Bq       -1.64444371      -0.05127407      10.78115585
  7-Bq       -1.08948640       2.22167614       9.91747987
 16-Bq       -1.97326041       1.60214544      11.18477703
  6-Bq        2.10261229      -2.06627375       7.93950154
  7-Bq        2.83337695      -2.96203946       7.81792923
  6-Bq        2.68182412       0.56874549       7.06480464
  7-Bq        3.67605422       0.86626543       6.54143620
  6-Bq       -2.15158210      -1.13761176      11.53471377
  7-Bq       -2.58370423      -2.01713951      12.15985160

@aug-cc-pVTZ.gbs/N

