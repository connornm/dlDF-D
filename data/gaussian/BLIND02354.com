%Chk=BLIND02354
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6         0.25261453        0.79397736        0.58897772
  6         0.73534221        0.89553906       -0.74611998
 16         0.00566299       -0.70971742        1.46443966
 16         1.13958875       -0.49570644       -1.76376876
  6        -0.45863173       -1.78495501        0.13311828
  6        -0.00617354       -1.69830486       -1.14571611
  6         0.07593644        2.05744678        1.11598161
  7         0.91966467        2.09933957       -1.23062960
 16         0.54455886        3.23185605       -0.06975970
  6        -1.34941218       -2.82187443        0.53610458
  7        -2.05571460       -3.66275457        0.91722204
  6        -0.40721996       -2.64591611       -2.13320617
  7        -0.68401956       -3.41197197       -2.96234219
  6        -0.37524803        2.41519633        2.40966483
  7        -0.73499131        2.72696348        3.47008075
  6        -0.25261414       -0.79397746        9.73113224
  6        -0.73534197       -0.89553934       11.06622986
 16        -0.00566304        0.70971742        8.85567034
 16        -1.13958921        0.49570600       12.08387859
  6         0.45863106        1.78495519       10.18699179
  6         0.00617271        1.69830487       11.46582611
  6        -0.07593548       -2.05744681        9.20412838
  7        -0.91966405       -2.09933992       11.55073946
 16        -0.54455762       -3.23185626       10.38986962
  6         1.34941118        2.82187495        9.78400563
  7         2.05571332        3.66275536        9.40288827
  6         0.40721861        2.64591627       12.45331623
  7         0.68401779        3.41197223       13.28245229
  6         0.37524933       -2.41519619        7.91044523
  7         0.73499289       -2.72696320        6.85002936

@blind-aug.gbs/N
