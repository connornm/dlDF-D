%Chk=BLIND01244
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.93703844       -0.20433850       -0.34830109
  6         0.31964428       -1.30385036        0.31179342
 16         0.65534797        1.48933536        0.02671354
 16        -0.87071845       -1.13950657        1.61203034
  6        -1.04520612        1.43023329        0.52545677
  6        -1.64710621        0.38676226        1.15516190
  6         1.84363140       -0.67113947       -1.27861785
  7         0.68729892       -2.50848318       -0.05090462
 16         1.86143943       -2.40369895       -1.22620510
  6        -1.75352531        2.63398653        0.24199296
  7        -2.28211132        3.64346613        0.01233933
  6        -3.01349870        0.45986219        1.55711949
  7        -4.11463458        0.48645771        1.92822402
  6         2.67995165        0.07622426       -2.14301610
  7         3.38005039        0.67365084       -2.85297326
  6         0.28470984       -0.40696572        6.82427657
  6        -1.11799506       -0.19135807        6.93271541
 16         1.52183614        0.56600285        7.60617638
 16        -1.85263416        1.08851870        7.91088164
  6         0.67777912        1.02269622        9.09718979
  6        -0.66039480        1.23047817        9.21449818
  6         0.52791127       -1.45137960        5.95516442
  7        -1.90402515       -0.97407463        6.23459544
 16        -0.98445880       -2.03880705        5.34484294
  6         1.55683954        1.20697031       10.20370934
  7         2.31530196        1.36844585       11.06965011
  6        -1.23452329        1.64283196       10.45315987
  7        -1.74116914        2.00715626       11.43372293
  6         1.77500775       -1.98815828        5.55293630
  7         2.78991610       -2.43516664        5.20474025

@blind-aug.gbs/N
