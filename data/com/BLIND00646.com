%Chk=BLIND00646
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.14330672        0.88860694        0.48057237
  6        -1.11015918        0.62191091       -0.52932673
 16         0.68050274       -0.33126685        1.44065703
 16        -1.62156140       -1.00379174       -1.00906314
  6         0.78103456       -1.65365954        0.26367232
  6        -0.13419991       -1.91854554       -0.70577357
  6        -0.00365869        2.25247246        0.64079870
  7        -1.67469719        1.65209508       -1.11059899
 16        -1.08934760        3.06000515       -0.44271278
  6         1.92565321       -2.48216537        0.44945683
  7         2.86009854       -3.14223170        0.65516875
  6         0.01816460       -3.03961419       -1.57406619
  7         0.08653383       -3.95741066       -2.28392632
  6         0.83789390        2.95032326        1.54077755
  7         1.51529623        3.53938580        2.27917264
  6        -0.61337976       -0.64581422        6.86313332
  6         0.72118678       -1.11949446        6.72035433
 16        -1.07939301        1.04371127        6.99294810
 16         2.14178624       -0.06790795        6.61714914
  6         0.16215896        1.79919889        5.97711876
  6         1.43909390        1.35671193        5.83121980
  6        -1.47868923       -1.71722285        6.95601115
  7         0.90107945       -2.41757434        6.69828006
 16        -0.56297466       -3.18734816        6.88512416
  6        -0.29002766        2.98603864        5.33056766
  7        -0.69700757        3.96025888        4.84430097
  6         2.37567388        2.06738647        5.02393945
  7         3.17960019        2.62964283        4.40064465
  6        -2.88535977       -1.69974505        7.11769285
  7        -4.03901658       -1.70427707        7.25964826

@aug-cc-pVTZ.gbs/N

