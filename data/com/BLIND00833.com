%Chk=BLIND00833
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.84253833        0.12346930       -0.56213291
  6        -0.61593783       -1.22260739       -0.15895008
 16        -0.36276633        1.54751967        0.34913875
 16         0.22479392       -1.68472819        1.32906969
  6         1.12202529        0.95422476        1.11563538
  6         1.35147729       -0.32793667        1.50429704
  6        -1.55911807        0.14024671       -1.74176641
  7        -1.09268051       -2.17073563       -0.92812450
 16        -1.90053052       -1.49186959       -2.21553903
  6         2.08297460        1.98168039        1.34367270
  7         2.82456108        2.85714521        1.53023062
  6         2.56588243       -0.69197593        2.15747199
  7         3.52704120       -1.02750185        2.71829488
  6        -2.00684391        1.26093300       -2.48255805
  7        -2.39009395        2.16921736       -3.09837852
  6         1.00633470       -0.04013000        7.44462802
  6         0.23113425       -0.06553119        8.63804568
 16         0.49690111       -0.68851973        5.89266642
 16        -1.40468890       -0.73431394        8.74481934
  6        -1.24861894       -0.38776672        5.97530707
  6        -2.00013111       -0.40795169        7.10776726
  6         2.24620934        0.50629774        7.70812531
  7         0.78536433        0.40828048        9.72716989
 16         2.33840609        0.90525394        9.39251253
  6        -1.83756016       -0.15654701        4.69828893
  7        -2.27570151        0.00272683        3.63346170
  6        -3.41000813       -0.19881270        7.05827376
  7        -4.56431045       -0.06279824        7.06178978
  6         3.32578800        0.70844423        6.81436072
  7         4.22535017        0.87317677        6.09690201

@aug-cc-pVTZ.gbs/N

