%Chk=BLIND01862
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.36929127        0.84887647       -0.42911725
  6        -0.23444949       -0.23904494       -1.33690102
 16        -0.70727735        0.69778712        1.28886240
 16        -0.35538798       -1.94241642       -0.86947424
  6         0.12503950       -0.82614716        1.64802122
  6         0.26174770       -1.87181950        0.79039410
  6        -0.26541745        2.04182659       -1.11569095
  7        -0.04193163        0.05052056       -2.60060538
 16        -0.04051552        1.70304919       -2.80076530
  6         0.62385668       -0.89039768        2.98149148
  7         0.99128906       -0.90537311        4.08419825
  6         0.91057024       -3.07532138        1.19604724
  7         1.40526780       -4.08446819        1.49239698
  6        -0.35164357        3.36004227       -0.60574888
  7        -0.42898855        4.44876738       -0.20600984
  6        -0.35963473       -0.81132806        7.15905374
  6        -0.57378609       -0.86939234        8.56496782
 16        -0.29194450        0.66333929        6.20548630
 16        -0.77719550        0.55438593        9.59747370
  6         0.41852550        1.77848069        7.38693207
  6         0.22310354        1.73345933        8.73133803
  6        -0.28465246       -2.09127448        6.64776319
  7        -0.65683718       -2.05679738        9.11363147
 16        -0.51061208       -3.22638969        7.93822135
  6         1.21094106        2.79949683        6.78635459
  7         1.82723994        3.62582156        6.24918038
  6         0.80520295        2.71042104        9.59194713
  7         1.23525355        3.50131824       10.32710985
  6        -0.09197882       -2.49089402        5.30301432
  7         0.05611580       -2.83694148        4.20328023

@blind-aug.gbs/N

