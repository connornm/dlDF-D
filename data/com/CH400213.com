%Chk=CH400213
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.95870203        0.00000000        0.55396424
  1         0.08290172        0.66531344       -0.88117608
  1        -0.23253422       -1.02932198       -0.33527607
  1        -0.80906953        0.36400854        0.66248792
  0        -0.63438720        0.00000000       -0.36656626
  0        -0.05485728       -0.44024767        0.58308714
  0         0.15387131        0.68111745        0.22185710
  0         0.53537317       -0.24086979       -0.43837797
  6         0.00000000        0.00000000        3.57785565
  1        -0.79503206       -0.51105275        3.00102373
  1         0.75338985       -0.74228592        3.90557703
  1         0.48671200        0.76278646        2.93969093
  1        -0.44506979        0.49055221        4.46513090
  0         0.52608438        0.33817110        3.95955379
  0        -0.49852911        0.49118148        3.36099760
  0        -0.32206447       -0.50474699        4.00013860
  0         0.29450921       -0.32460559        2.99073260

@aug-cc-pVTZ.gbs/N

