%Chk=CH400270
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.08203946        0.00000000        0.23489693
  1        -0.53016002        0.67197299        0.70240202
  1        -0.15248234        0.35586840       -1.03734892
  1        -0.39939709       -1.02784139        0.10004997
  0        -0.71600138        0.00000000       -0.15543474
  0         0.35081466       -0.44465439       -0.46478972
  0         0.10089980       -0.23548334        0.68642899
  0         0.26428692        0.68013772       -0.06620453
  6         0.00000000        0.00000000        2.91107902
  1         0.95020398        0.50654429        2.65318129
  1        -0.77169498        0.75921198        3.14360855
  1        -0.33452422       -0.61572403        2.05377895
  1         0.15601522       -0.65003224        3.79374731
  0        -0.62876391       -0.33518778        3.08173374
  0         0.51064189       -0.50238171        2.75721082
  0         0.22135958        0.40743362        3.47836706
  0        -0.10323756        0.43013587        2.32700447

@aug-cc-pVTZ.gbs/N

