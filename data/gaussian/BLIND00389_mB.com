%Chk=BLIND00389_mB
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.20709165       -0.43412408        6.75273306
  6        -0.05781367       -1.36710556        5.68826686
 16         0.61890127        1.11293754        6.86610261
 16         0.97276320       -1.10180852        4.27330498
  6         0.73577949        1.54215125        5.14964197
  6         0.87739179        0.66159088        4.12387228
  6        -1.04351112       -0.95982028        7.71676065
  7        -0.70222831       -2.50402848        5.78824734
 16        -1.53529886       -2.54877267        7.22870155
  6         0.71348879        2.94882852        4.92193196
  7         0.71491563        4.10356141        4.78857928
  6         1.01029520        1.11548609        2.77841409
  7         1.15111106        1.44098382        1.67154993
  6        -1.45101946       -0.37555133        8.94058057
  7        -1.78832151        0.08480483        9.95319869

@blind-aug.gbs/N

