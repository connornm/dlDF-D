%Chk=CH400224
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.54684026       0.00000000       0.96278330
  1-Bq        0.54810279       0.61984871      -0.73576957
  1-Bq       -0.08070238      -1.03736140      -0.37862172
  1-Bq       -1.01424067       0.41751269       0.15160799
  6         0.00000000        0.00000000        4.12585674
  1        -0.26900754       -0.80851210        3.41880010
  1        -0.84307486        0.17781010        4.82126619
  1         0.89746764       -0.29720398        4.70223126
  1         0.21461477        0.92790597        3.56112941

@aug-cc-pVTZ.gbs/N

