%Chk=CH400260
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.68175711        0.00000000        0.87246387
  1        -1.03343542        0.20731437        0.33914294
  1         0.31730928        0.78239560       -0.71635041
  1         0.03436903       -0.98970997       -0.49525640
  0        -0.45112868        0.00000000       -0.57732214
  0         0.68383938       -0.13718296       -0.22441586
  0        -0.20996821       -0.51772264        0.47401957
  0        -0.02274249        0.65490560        0.32771843
  6         0.00000000        0.00000000        3.71317576
  1        -0.32777287       -0.52723081        2.79634504
  1        -0.10629777       -0.67451231        4.58479487
  1         1.06015601        0.30119216        3.60669354
  1        -0.62608538        0.90055096        3.86486960
  0         0.21689211        0.34887636        4.31985610
  0         0.07033879        0.44633469        3.13641261
  0        -0.70152079       -0.19930327        3.78363661
  0         0.41428988       -0.59590778        3.61279772

@aug-cc-pVTZ.gbs/N

