%Chk=BLIND00645
%Mem=1GB
%NProcShared=48
#T M05/Gen test Massage SCF=(fermi,maxcyc=60)

M05 opt

0 1
  6        -0.63293055        0.77764251        0.18915588
  6         0.11334924        0.92692629       -1.01359409
 16        -1.13954740       -0.75575467        0.88237237
 16         0.68612882       -0.42670718       -2.00060387
  6         0.21201807       -1.79328641        0.39153559
  6         0.93319321       -1.66102997       -0.75304237
  6        -0.95958006        2.02169075        0.69010228
  7         0.35749747        2.14748607       -1.42406394
 16        -0.34104102        3.23831424       -0.37857134
  6         0.46675041       -2.85235069        1.31054808
  7         0.62075061       -3.71273466        2.07684853
  6         1.97352593       -2.58122038       -1.07692408
  7         2.80934421       -3.32360087       -1.39505163
  6        -1.70119539        2.33312493        1.85551375
  7        -2.31567465        2.60702996        2.80339940
  6        -0.82106076        0.59575072        5.36406687
  6        -1.13432687       -0.77339260        5.13355100
 16         0.61557040        1.18020153        6.19052001
 16        -0.09830192       -2.12406257        5.62030691
  6         1.78780588       -0.06301004        5.71679959
  6         1.50181366       -1.37283084        5.49269702
  6        -1.85081256        1.39134286        4.90368921
  7        -2.27898256       -1.04061114        4.55370825
 16        -3.10184088        0.37687106        4.26290715
  6         3.12184236        0.43058122        5.62845931
  7         4.19580127        0.87390168        5.59279066
  6         2.53026454       -2.30354792        5.16103991
  7         3.33926186       -3.10173387        4.91742409
  6        -1.95645124        2.80287760        4.94342065
  7        -2.06413425        3.95978647        4.97601141

@blind-aug.gbs/N

