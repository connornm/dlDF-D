%Chk=BLIND01971
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.74668729       -0.21936192       -0.65988434
  6         0.40223940       -1.30361901        0.19543364
 16         0.63037779        1.48185447       -0.23461791
 16        -0.23716772       -1.10985068        1.83500697
  6        -0.77565984        1.44847175        0.84522739
  6        -1.11627174        0.41970227        1.66575218
  6         1.25311152       -0.70780344       -1.84742912
  7         0.60631714       -2.51671057       -0.25691709
 16         1.27745262       -2.43969641       -1.77828584
  6        -1.53075370        2.65685284        0.81806034
  7        -2.10006997        3.66955236        0.77917992
  6        -2.24480182        0.51346547        2.53268370
  7        -3.13742293        0.55728541        3.27580785
  6         1.72565444        0.01935815       -2.96685204
  7         2.12597433        0.60004135       -3.89078218
  6         0.74668718        0.21936192        7.64732553
  6         0.40223944        1.30361901        8.50264357
 16         0.63037774       -1.48185447        8.07259198
 16        -0.23716740        1.10985068       10.14221701
  6        -0.77565970       -1.44847175        9.15243752
  6        -1.11627145       -0.41970226        9.97296238
  6         1.25311120        0.70780344        6.45978066
  7         0.60631710        2.51671057        8.05029280
 16         1.27745232        2.43969641        6.52892394
  6        -1.53075356       -2.65685284        9.12527061
  7        -2.10006984       -3.66955236        9.08639028
  6        -2.24480138       -0.51346547       10.83989409
  7        -3.13742236       -0.55728540       11.58301839
  6         1.72565392       -0.01935815        5.34035766
  7         2.12597366       -0.60004136        4.41642745

@aug-cc-pVTZ.gbs/N
