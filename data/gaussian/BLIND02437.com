%Chk=BLIND02437
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.47484618       -0.62657373       -0.65041170
  6         0.40508425       -1.12844336        0.67965852
 16         0.55931000        1.07346591       -1.08773051
 16         0.35937994       -0.10663717        2.12479227
  6        -0.40291225        1.79005472        0.21793059
  6        -0.47907405        1.32046294        1.49124398
  6         0.53824200       -1.67929789       -1.54094863
  7         0.40884774       -2.43000193        0.83390233
 16         0.53225660       -3.16805287       -0.65303007
  6        -1.08657501        2.97644635       -0.17725264
  7        -1.60515319        3.95162574       -0.53959898
  6        -1.24768443        2.00118063        2.48109419
  7        -1.83726038        2.53870914        3.32630980
  6         0.62918926       -1.63181075       -2.95326659
  7         0.71336500       -1.61159504       -4.11240493
  6         0.35958182        0.73381438       10.21734241
  6         0.35721666       -0.55802400       10.81481599
 16         0.56033894        1.05502274        8.50123099
 16         0.52633609       -2.07029833        9.90982246
  6        -0.21941051       -0.38040289        7.81154545
  6        -0.22999636       -1.61874130        8.37194706
  6         0.25203751        1.69900963       11.19822750
  7         0.25617605       -0.61835475       12.12019372
 16         0.18631615        0.92233419       12.74654353
  6        -0.82380459       -0.13441841        6.54454541
  7        -1.28072680        0.10898630        5.50379124
  6        -0.84752899       -2.72030198        7.70925598
  7        -1.31257253       -3.65172515        7.19242767
  6         0.22735860        3.10677907       11.04741437
  7         0.21486491        4.26438530       10.94307050

@blind-aug.gbs/N

