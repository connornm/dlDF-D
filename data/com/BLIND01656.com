%Chk=BLIND01656
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.68352633       -0.64116097       -0.40350127
  6        -0.09954238       -1.11399532        0.80530685
 16        -0.94705192        1.04888720       -0.80699609
 16         0.49211949       -0.06105201        2.09984730
  6         0.43859480        1.79489873        0.01005473
  6         1.00619969        1.35311333        1.16343862
  6        -1.07653483       -1.71302059       -1.17936498
  7        -0.03077069       -2.41190964        0.97427101
 16        -0.71298872       -3.18217516       -0.33423443
  6         0.90453995        2.97336899       -0.64201321
  7         1.23302190        3.94116574       -1.19578400
  6         2.09144121        2.05615738        1.76503802
  7         2.95700878        2.61267030        2.30541964
  6        -1.70750987       -1.69627252       -2.44694474
  7        -2.23398092       -1.70129388       -3.48323569
  6         0.07250503        0.65660336        7.04140001
  6         0.80642734        1.05255291        8.19483343
 16        -0.24043982       -1.00488529        6.56177765
 16         1.50383480       -0.07918195        9.36418684
  6        -0.35030487       -1.78359454        8.15101132
  6         0.34373335       -1.41497134        9.26009740
  6        -0.29465165        1.77388657        6.31878156
  7         1.00023098        2.33626266        8.37480408
 16         0.32045307        3.18554563        7.11482894
  6        -1.22947209       -2.90518379        8.16487552
  7        -1.93768893       -3.82587798        8.12078334
  6         0.21736389       -2.14060637       10.48133714
  7         0.16774550       -2.71942391       11.48803193
  6        -1.02146656        1.83619068        5.10509843
  7        -1.60866652        1.90561412        4.10436245

@aug-cc-pVTZ.gbs/N

