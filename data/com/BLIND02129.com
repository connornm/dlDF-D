%Chk=BLIND02129
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.63319889       -0.51375040       -0.61337477
  6         0.29502387       -1.21608591        0.57750270
 16         0.73692415        1.23422359       -0.76285602
 16        -0.11191415       -0.43075604        2.11120631
  6        -0.51917926        1.73078738        0.38594801
  6        -0.85253914        1.06906229        1.52552475
  6         0.93688545       -1.41529306       -1.61352401
  7         0.32197412       -2.52573594        0.53333816
 16         0.80136064       -3.02352533       -0.98103008
  6        -1.15074680        2.95740664        0.02871239
  7        -1.61936398        3.97187659       -0.29135273
  6        -1.84899025        1.58093439        2.40814249
  7        -2.63289125        1.97546035        3.17022819
  6         1.33373691       -1.14890889       -2.94645943
  7         1.66956047       -0.94871805       -4.04110188
  6        -0.64809836       -0.34736206       11.51842840
  6        -0.67364068       -1.09059916       12.73200693
 16        -0.28230585        1.36666010       11.38865740
 16        -0.30505996       -0.40631631       14.32275945
  6         0.88049008        1.54911005       12.71483225
  6         0.86800511        0.84542914       13.87772577
  6        -1.02392816       -1.16046657       10.46823021
  7        -1.02542886       -2.35110640       12.65971168
 16        -1.39311753       -2.73857213       11.08313557
  6         1.84714023        2.56851289       12.47505109
  7         2.59874923        3.42478857       12.24468422
  6         1.82448754        1.10370993       14.90357874
  7         2.56515946        1.30024657       15.77748067
  6        -1.14476821       -0.82653921        9.09744261
  7        -1.25820571       -0.56775139        7.96994240

@aug-cc-pVTZ.gbs/N
