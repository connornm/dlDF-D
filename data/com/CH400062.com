%Chk=CH400062
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.85381971        0.00000000       -0.70496657
  1        -0.86715396       -0.50259373       -0.47056297
  1         0.28379070       -0.54108791        0.92340277
  1        -0.27045645        1.04368164        0.25212677
  0        -0.56498503        0.00000000        0.46648672
  0         0.57380850        0.33257364        0.31137842
  0        -0.18778847        0.35804581       -0.61102915
  0         0.17896500       -0.69061945       -0.16683598
  6         0.00000000        0.00000000        3.27633452
  1        -0.83267198        0.37763607        2.65180723
  1         0.28239893        0.77026223        4.01992910
  1        -0.32074767       -0.92028698        3.80186164
  1         0.87102072       -0.22761132        2.63174011
  0         0.55099126       -0.24988733        3.68959339
  0        -0.18686752       -0.50969382        2.78428708
  0         0.21224344        0.60896740        2.92858551
  0        -0.57636719        0.15061375        3.70287210

@aug-cc-pVTZ.gbs/N

