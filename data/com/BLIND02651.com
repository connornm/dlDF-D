%Chk=BLIND02651
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.69415776       -0.56349400       -0.49165884
  6         0.18367347       -1.18224402        0.68408810
 16         0.86378892        1.16918860       -0.73156075
 16        -0.39979402       -0.29351831        2.09982842
  6        -0.51842920        1.76125940        0.20814523
  6        -1.01630786        1.17967941        1.33132180
  6         1.10209773       -1.53232564       -1.38622771
  7         0.18080754       -2.49232770        0.72323594
 16         0.84129645       -3.09427617       -0.68116049
  6        -1.06425003        2.97433164       -0.30334770
  7        -1.45920320        3.97502463       -0.74359599
  6        -2.10584742        1.76523160        2.04123315
  7        -2.97209640        2.22204450        2.66722101
  6         1.67760981       -1.35869927       -2.66833854
  7         2.15956313       -1.23482900       -3.71880181
  6        -0.69415778        0.56349405       11.33858943
  6        -0.18367319        1.18224410       10.16284264
 16        -0.86378920       -1.16918855       11.57849114
 16         0.39979444        0.29351846        8.74710234
  6         0.51842901       -1.76125944       10.63878536
  6         1.01630795       -1.17967942        9.51560893
  6        -1.10209779        1.53232565       12.23315832
  7        -0.18080708        2.49232779       10.12369492
 16        -0.84129617        3.09427621       11.52809129
  6         1.06424958       -2.97433181       11.15027827
  7         1.45920254       -3.97502489       11.59052654
  6         2.10584756       -1.76523169        8.80569773
  7         2.97209660       -2.22204464        8.17970998
  6        -1.67761012        1.35869924       13.51526903
  7        -2.15956365        1.23482893       14.56573219

@aug-cc-pVTZ.gbs/N

