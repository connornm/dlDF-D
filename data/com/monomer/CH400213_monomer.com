%Chk=CH400213
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.95870203        0.00000000        0.55396424
  1         0.08290172        0.66531344       -0.88117608
  1        -0.23253422       -1.02932198       -0.33527607
  1        -0.80906953        0.36400854        0.66248792
  6-Bq        0.00000000       0.00000000       3.57785565
  1-Bq       -0.79503206      -0.51105275       3.00102373
  1-Bq        0.75338985      -0.74228592       3.90557703
  1-Bq        0.48671200       0.76278646       2.93969093
  1-Bq       -0.44506979       0.49055221       4.46513090

@aug-cc-pVTZ.gbs/N

