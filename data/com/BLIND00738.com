%Chk=BLIND00738
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.39904215       -0.41411546        0.84284202
  6         0.39123195        1.00914065        0.85322101
 16        -0.78199933       -1.42693213        0.02561497
 16        -0.80652435        2.00118336        0.00715988
  6        -1.17313180       -0.40425301       -1.36911256
  6        -1.18326738        0.95498466       -1.37288948
  6         1.43293542       -0.87831325        1.63080669
  7         1.31309760        1.61193590        1.56361641
 16         2.26672950        0.47628041        2.31978418
  6        -1.53601139       -1.15194679       -2.52685263
  7        -1.84760590       -1.80408990       -3.43724335
  6        -1.55977440        1.68410454       -2.53932362
  7        -1.88744992        2.31822675       -3.45662614
  6         1.79435220       -2.21939209        1.90658540
  7         2.09945202       -3.31400138        2.15120861
  6        -0.99544019        0.18627243        6.30552119
  6        -0.51331232        0.49939526        5.00347147
 16        -0.35881870       -1.08141471        7.34292132
 16         0.84345987       -0.32729215        4.22222535
  6         1.35089555       -1.05328404        6.87365814
  6         1.82442955       -0.75511742        5.63490624
  6        -2.08445767        0.98135032        6.60080970
  7        -1.14140093        1.43733971        4.33741973
 16        -2.41644615        2.00697001        5.24334351
  6         2.22794242       -1.42802222        7.93264337
  7         2.90074163       -1.76257103        8.81955449
  6         3.22070484       -0.80784862        5.34928996
  7         4.34597265       -0.86608929        5.06414698
  6        -2.87136493        0.99257470        7.77801801
  7        -3.53358250        1.00905356        8.73315793

@aug-cc-pVTZ.gbs/N
