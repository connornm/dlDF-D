%Chk=CH400220
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.06486086        0.00000000        0.30340946
  1        -0.48155146       -0.93612142        0.34317755
  1        -0.51380617        0.86816255        0.45638028
  1        -0.06950324        0.06795887       -1.10296729
  6-Bq        0.00000000       0.00000000       2.69548980
  1-Bq        0.59062195      -0.76236902       2.15148139
  1-Bq        0.39167093       1.00827671       2.45893555
  1-Bq       -1.06121569      -0.06557706       2.38644794
  1-Bq        0.07892281      -0.18033064       3.78509432

@aug-cc-pVTZ.gbs/N

