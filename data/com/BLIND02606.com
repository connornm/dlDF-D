%Chk=BLIND02606
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.59513136       -0.74088939        0.37149763
  6         0.82456742        0.53065563        0.96845576
 16         0.05872949       -1.00103457       -1.28171566
 16         0.57261740        2.07545204        0.14098154
  6        -0.96154453        0.42978756       -1.51817085
  6        -0.75434974        1.64879675       -0.95360823
  6         0.92630662       -1.73959107        1.26492984
  7         1.28232280        0.54637942        2.19649191
 16         1.50724495       -1.01480584        2.72849521
  6        -2.03979179        0.20479049       -2.42252078
  7        -2.88899377       -0.02006599       -3.18376226
  6        -1.61092242        2.75043365       -1.24792325
  7        -2.26507441        3.68237348       -1.48137441
  6         0.86247264       -3.14195066        1.07923859
  7         0.82534591       -4.29548104        0.94108940
  6        -0.59513134        0.74088939        7.10263028
  6        -0.82456742       -0.53065565        6.50567219
 16        -0.05872949        1.00103462        8.75584356
 16        -0.57261743       -2.07545203        7.33314647
  6         0.96154451       -0.42978751        8.99229881
  6         0.75434970       -1.64879673        8.42773624
  6        -0.92630658        1.73959104        6.20919802
  7        -1.28232278       -0.54637948        5.27763603
 16        -1.50724490        1.01480576        4.74563267
  6         2.03979176       -0.20479043        9.89664875
  7         2.88899374        0.02006606       10.65789022
  6         1.61092236       -2.75043363        8.72205130
  7         2.26507434       -3.68237346        8.95550251
  6        -0.86247258        3.14195064        6.39488921
  7        -0.82534583        4.29548101        6.53303837

@aug-cc-pVTZ.gbs/N

