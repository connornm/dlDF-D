%Chk=BLIND02178
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.24949176        0.15745258       -0.97676581
  6         0.34182663       -1.20489810       -0.57514216
 16         0.51434653        1.54224594        0.07238429
 16         0.71172824       -1.73028743        1.07462363
  6        -0.07824004        0.89647778        1.61379491
  6         0.00270671       -0.40159297        2.00889489
  6        -0.00908979        0.22558717       -2.33094323
  7         0.17337343       -2.11869281       -1.49951371
 16        -0.08493712       -1.38374109       -2.97058085
  6        -0.63123365        1.89755083        2.46417158
  7        -1.05030995        2.75226081        3.13128396
  6        -0.46349719       -0.80947212        3.29350785
  7        -0.80343035       -1.18106287        4.34102487
  6        -0.16517980        1.37771228       -3.13928232
  7        -0.28690026        2.31235951       -3.81951361
  6        -0.18657700        0.18053491       11.11697495
  6         0.96303848       -0.62752569       10.89059264
 16        -0.78029017        1.42433126       10.02656588
 16         1.97446825       -0.53040514        9.44070790
  6        -0.32293443        0.72615795        8.46211650
  6         0.77087841       -0.04753411        8.23273012
  6        -0.72263530       -0.11146847       12.35498899
  7         1.30390866       -1.46926957       11.83565073
 16         0.25314553       -1.32380432       13.11848103
  6        -1.20710235        1.09241107        7.40608837
  7        -1.94197264        1.43732303        6.57410616
  6         1.07454649       -0.52124935        6.92226397
  7         1.37641657       -0.90451451        5.86731924
  6        -1.85579619        0.46329889       12.98009484
  7        -2.77771365        0.93045012       13.51201010

@aug-cc-pVTZ.gbs/N

