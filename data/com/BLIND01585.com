%Chk=BLIND01585
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.53234182        0.38284131       -0.78176315
  6        -0.53875327       -1.03941234       -0.72717281
 16        -0.41962538        1.44879703        0.61091415
 16        -0.40029048       -1.97575213        0.76905050
  6         0.59740569        0.47003592        1.68413201
  6         0.60179996       -0.88785526        1.74542051
  6        -0.70000552        0.79793219       -2.08747815
  7        -0.69038398       -1.68545367       -1.85743771
 16        -0.87369264       -0.59625399       -3.10281270
  6         1.41525415        1.25420638        2.54859151
  7         2.04495833        1.93577468        3.24868459
  6         1.42709103       -1.57885350        2.68102621
  7         2.06059986       -2.18200642        3.44640144
  6        -0.76427066        2.12100703       -2.58798928
  7        -0.82858533        3.20017754       -3.01501452
  6         0.38047709       -0.29624338        8.60141593
  6         0.41472304        1.11259906        8.80095756
 16         0.50256634       -1.50210744        9.87401276
 16         0.55553340        1.88716856       10.38677056
  6        -0.29561128       -0.64635464       11.20619958
  6        -0.27163934        0.69771645       11.40755075
  6         0.30884513       -0.57208521        7.25082962
  7         0.37202104        1.87372946        7.73478822
 16         0.31781652        0.92126228        6.37068227
  6        -0.95739469       -1.52101675       12.11619879
  7        -1.46182382       -2.27542865       12.84256437
  6        -0.90919824        1.28271849       12.54123201
  7        -1.38930746        1.79919439       13.46517676
  6         0.26643389       -1.83543066        6.61260711
  7         0.24022140       -2.86385755        6.07154615

@aug-cc-pVTZ.gbs/N

