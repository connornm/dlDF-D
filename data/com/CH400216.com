%Chk=CH400216
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.98511389        0.00000000       -0.50550627
  1         0.03758977       -0.66874685        0.88167489
  1        -0.24694286        1.02856879        0.32718713
  1        -0.77576079       -0.35982195       -0.70335576
  0        -0.65186431        0.00000000        0.33450091
  0        -0.02487370        0.44251960       -0.58341721
  0         0.16340571       -0.68061906       -0.21650452
  0         0.51333229        0.23809946        0.46542082
  6         0.00000000        0.00000000        2.85830957
  1        -0.97902404        0.13732203        3.35694826
  1         0.73696195       -0.38222146        3.59096214
  1        -0.10672971       -0.72599283        2.02913580
  1         0.34879179        0.97089225        2.45619208
  0         0.64783457       -0.09086800        2.52835303
  0        -0.48765854        0.25292155        2.37350262
  0         0.07062461        0.48040010        3.40698601
  0        -0.23080064       -0.64245365        3.12439661

@aug-cc-pVTZ.gbs/N

