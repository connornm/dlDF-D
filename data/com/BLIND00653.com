%Chk=BLIND00653
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.42666599        0.33646733        0.86362884
  6        -0.42627714        1.28881172        0.23804524
 16         0.22644073       -1.40771149        0.78446628
 16        -1.82559303        0.86333649       -0.75968930
  6        -0.46762292       -1.57887596       -0.83821683
  6        -1.28213561       -0.67707447       -1.44728759
  6         1.39253433        0.99187794        1.60040458
  7        -0.17002307        2.55699712        0.44756771
 16         1.13630678        2.70124735        1.46914056
  6        -0.12582290       -2.81198006       -1.46568530
  7         0.16037034       -3.84017036       -1.92624058
  6        -1.82663774       -0.93397400       -2.74013793
  7        -2.31581716       -1.11412205       -3.77897203
  6         2.43009925        0.43235902        2.38498476
  7         3.28184248       -0.00960212        3.04095712
  6        -0.78382294       -0.52170147        5.92384446
  6         0.44816616       -1.23398144        5.95005076
 16        -0.94591507        1.22473745        6.03427023
 16         2.03776901       -0.46329229        6.06747879
  6         0.54073489        1.73854177        5.21573550
  6         1.72276595        1.06758319        5.23189719
  6        -1.83510694       -1.41460318        5.87063057
  7         0.38860426       -2.54290866        5.91808307
 16        -1.20319836       -3.02870008        5.88509797
  6         0.40685938        2.98924401        4.54572736
  7         0.25600143        4.02246280        4.03494967
  6         2.87565324        1.59355406        4.57756740
  7         3.84692204        1.99801795        4.08363216
  6        -3.22342527       -1.13732727        5.84132638
  7        -4.36677270       -0.92847452        5.82605919

@aug-cc-pVTZ.gbs/N
