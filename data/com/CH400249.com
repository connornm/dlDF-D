%Chk=CH400249
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.99499086        0.00000000       -0.48577684
  1         0.11879143        0.18864048        1.08456873
  1        -0.62856480        0.79485638       -0.44620128
  1        -0.48521749       -0.98349685       -0.15259060
  0        -0.65840005        0.00000000        0.32144566
  0        -0.07860603       -0.12482617       -0.71767504
  0         0.41593055       -0.52596812        0.29525793
  0         0.32107553        0.65079429        0.10097144
  6         0.00000000        0.00000000        4.67338525
  1         0.02257061        0.06184028        3.56810144
  1         0.16911481       -1.04806721        4.98793454
  1         0.79499034        0.64486500        5.09543435
  1        -0.98667576        0.34136192        5.04207067
  0        -0.01493531       -0.04092062        5.40476776
  0        -0.11190575        0.69352145        4.46524337
  0        -0.52605677       -0.42671663        4.39410917
  0         0.65289782       -0.22588419        4.42942070

@aug-cc-pVTZ.gbs/N

