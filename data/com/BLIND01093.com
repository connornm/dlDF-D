%Chk=BLIND01093
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.11174045        0.84749992        0.55710580
  6         1.00433008        0.76790517       -0.54868581
 16        -0.52124824       -0.52869294        1.44813619
 16         1.60309172       -0.75105334       -1.23340320
  6        -0.61637207       -1.73523337        0.15224509
  6         0.22839636       -1.82160141       -0.90914593
  6        -0.12467876        2.17133335        0.86835201
  7         1.42892636        1.89612550       -1.06317285
 16         0.79122996        3.17438927       -0.20859581
  6        -1.66826743       -2.67795459        0.34167185
  7        -2.52358563       -3.43661179        0.55158389
  6         0.09253211       -2.86300374       -1.87401070
  7         0.03734826       -3.71057324       -2.66742536
  6        -0.93796259        2.70022704        1.89986028
  7        -1.59362121        3.15170418        2.74684085
  6        -0.30329094       -0.44491328        6.48404476
  6        -1.26139188       -0.52427911        7.53360196
 16         0.85103109        0.86183066        6.26386856
 16        -1.42524821        0.68012188        8.82090833
  6         1.15607378        1.30994172        7.95214209
  6         0.25088174        1.23799639        8.96362382
  6        -0.46599894       -1.51330717        5.62534584
  7        -2.09058981       -1.53918382        7.51750107
 16        -1.79027205       -2.48432441        6.18055104
  6         2.47055775        1.81663661        8.16773856
  7         3.54427793        2.24506218        8.28926855
  6         0.58498357        1.66952168       10.28123675
  7         0.80388747        2.03660736       11.36209535
  6         0.26340641       -1.81539263        4.44981165
  7         0.84541964       -2.07370461        3.47737582

@aug-cc-pVTZ.gbs/N
