%Chk=CH400256
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.81746264        0.00000000        0.74682045
  1        -0.05262665        0.99172091       -0.48959783
  1         0.19686870       -0.77815675       -0.76269308
  1        -0.96170469       -0.21356415        0.50547046
  0        -0.54092702        0.00000000       -0.49418204
  0         0.03482382       -0.65623627        0.32397406
  0        -0.13027090        0.51491774        0.50468520
  0         0.63637410        0.14131853       -0.33447722
  6         0.00000000        0.00000000        4.85382173
  1         0.75552769        0.15853776        4.06007767
  1         0.46985807        0.14792125        5.84545607
  1        -0.82748826        0.72445907        4.72572730
  1        -0.39789750       -1.03091808        4.78402590
  0        -0.49994376       -0.10490676        5.37905382
  0        -0.31091197       -0.09788166        4.19764275
  0         0.54756112       -0.47938519        4.93858370
  0         0.26329461        0.68217361        4.90000666

@aug-cc-pVTZ.gbs/N

