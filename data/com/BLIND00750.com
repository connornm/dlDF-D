%Chk=BLIND00750
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.14749782       -0.41411545        0.92079406
  6         0.13709740        1.00914065        0.92857570
 16        -0.75797867       -1.42693214       -0.19402940
 16        -0.77636630        2.00118336       -0.21860505
  6        -0.74359198       -0.40425302       -1.64249152
  6        -0.75226750        0.95498465       -1.64895144
  6         0.91987443       -0.87831324        1.96638384
  7         0.82359942        1.61193591        1.86836976
 16         1.52780447        0.47628043        2.86099192
  6        -0.76833338       -1.15194680       -2.85551729
  7        -0.81298595       -1.80408991       -3.81671891
  6        -0.78766236        1.68410453       -2.87413442
  7        -0.84582239        2.31822673       -3.84646799
  6         1.18978066       -2.21939208        2.33220706
  7         1.41432565       -3.31400137        2.65237259
  6         0.79798095        0.16387044        7.10503107
  6         1.02449250       -0.84872147        6.13077403
 16        -0.25216785        1.55681632        6.89194437
 16         0.25016205       -0.88538883        4.53902775
  6        -1.50060199        0.86134701        5.84225075
  6        -1.29797265       -0.10612190        4.90919699
  6         1.58293743       -0.07595760        8.21481891
  7         1.88702071       -1.78962148        6.42845222
 16         2.53217152       -1.49919538        7.93522028
  6        -2.78313515        1.45373263        6.03007417
  7        -3.80345153        1.98334249        6.20229758
  6        -2.36507242       -0.56386976        4.08112277
  7        -3.19742421       -0.94593331        3.36548036
  6         1.67802621        0.68473163        9.40539796
  7         1.77565619        1.30195776       10.38549806

@aug-cc-pVTZ.gbs/N
