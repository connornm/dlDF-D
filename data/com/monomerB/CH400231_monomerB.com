%Chk=CH400231
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        1.10522421       0.00000000       0.06682367
  1-Bq       -0.38857034       0.98901727       0.31119751
  1-Bq       -0.30663509      -0.20518597      -1.04396341
  1-Bq       -0.41001878      -0.78383130       0.66594223
  6         0.00000000        0.00000000        4.93581916
  1         0.01273771       -0.46498099        3.93102247
  1         0.72204917        0.83896838        4.96344155
  1        -1.01669106        0.38073836        5.15345082
  1         0.28190418       -0.75472576        5.69536180

@aug-cc-pVTZ.gbs/N

