%Chk=CH400161
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.90628788        0.00000000        0.63610394
  1        -0.62885962       -0.87535622        0.25352126
  1        -0.57422680        0.93026068        0.17568325
  1         0.29679854       -0.05490447       -1.06530845
  0        -0.59970399        0.00000000       -0.42091931
  0         0.41612564        0.57923605       -0.16775873
  0         0.37997430       -0.61556714       -0.11625218
  0        -0.19639594        0.03633109        0.70493022
  6         0.00000000        0.00000000        3.18769560
  1        -0.26068006        0.09391244        2.11568240
  1         0.42397507        0.95668381        3.54961712
  1        -0.91038340       -0.24359948        3.76894195
  1         0.74708839       -0.80699676        3.31654094
  0         0.17249582       -0.06214324        3.89706247
  0        -0.28055052       -0.63305171        2.94820682
  0         0.60241405        0.16119335        2.80307637
  0        -0.49435934        0.53400160        3.10243675

@aug-cc-pVTZ.gbs/N

