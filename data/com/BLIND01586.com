%Chk=BLIND01586
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.26202713       -0.07037672       -0.98361467
  6        -0.40578556       -1.29025729       -0.26455169
 16        -0.44477364        1.53698248       -0.29694388
 16        -0.76495262       -1.38546274        1.46649519
  6         0.14997137        1.25344279        1.34926562
  6         0.02027006        0.09377463        2.04639643
  6        -0.02810050       -0.34214605       -2.31655239
  7        -0.29619108       -2.40641304       -0.94281816
 16        -0.03570408       -2.05978641       -2.54996265
  6         0.76331191        2.40226479        1.92815338
  7         1.23278403        3.37146522        2.36569136
  6         0.49400920       -0.01488619        3.38718202
  7         0.83842668       -0.13939066        4.49026922
  6         0.16192566        0.57311370       -3.38020558
  7         0.31068799        1.30995434       -4.26678978
  6         0.58918154       -0.66674189       10.26789104
  6         0.68546458        0.66438653        9.77326657
 16         0.25006335       -1.10484018       11.93569960
 16         0.44444283        2.11133396       10.76464944
  6        -0.81195836        0.23332635       12.41061254
  6        -0.73148133        1.50798096       11.94539954
  6         0.87308713       -1.56003926        9.25461880
  7         1.00553478        0.81702013        8.51147495
 16         1.25270896       -0.67402105        7.81389483
  6        -1.77222181       -0.13728095       13.39632750
  7        -2.52068043       -0.48086935       14.21666264
  6        -1.60747248        2.52394249       12.42955201
  7        -2.28025489        3.38915468       12.81646003
  6         0.90333039       -2.97457071        9.31253745
  7         0.94188912       -4.13586763        9.34413824

@aug-cc-pVTZ.gbs/N

