%Chk=BLIND01661
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.75234177        0.68897338       -0.02016080
  6         0.38436373        1.06034927       -0.79204629
 16        -1.28625077       -0.96066886        0.26643700
 16         1.47566519       -0.09578680       -1.57114760
  6         0.27282829       -1.80508669        0.28522745
  6         1.36582467       -1.46042585       -0.44567054
  6        -1.41223292        1.82272303        0.40932565
  7         0.60863668        2.34107531       -0.95728091
 16        -0.58549191        3.21925209       -0.19979866
  6         0.27767769       -2.95210309        1.13105220
  7         0.22571253       -3.89105456        1.81432579
  6         2.56043054       -2.23764946       -0.39292378
  7         3.54366121       -2.85738675       -0.40337690
  6        -2.59053398        1.91164700        1.18962532
  7        -3.56186423        2.00335855        1.82145003
  6         0.75398434        0.68717489        4.32017923
  6        -0.38241560        1.06191849        5.09088609
 16         1.28413612       -0.96398238        4.03533646
 16        -1.47709858       -0.09094845        5.87008783
  6        -0.27695861       -1.80465632        4.01603271
  6        -1.36968246       -1.45676646        4.74580753
  6         1.41693579        1.81898032        3.89027801
  7        -0.60372817        2.34331564        5.25490474
 16         0.59309509        3.21799407        4.49762765
  6        -0.28392245       -2.95234778        3.17113894
  7        -0.23369536       -3.89197917        2.48867055
  6        -2.56611713       -2.23115229        4.69277354
  7        -3.55084629       -2.84851002        4.70297395
  6         2.59604698        1.90442800        3.11081413
  7         3.56808064        1.99328314        2.47966339

@aug-cc-pVTZ.gbs/N
