%Chk=BLIND00199
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.33620027       -0.92324388       -0.27513454
  6        -1.19954692       -0.24692075        0.63208441
 16         0.69457952       -0.14182914       -1.46484067
 16        -1.35208291        1.51419796        0.73098835
  6         1.09721728        1.35701636       -0.60724688
  6         0.28288803        2.01204241        0.26192099
  6        -0.48470701       -2.28759207       -0.12747617
  7        -1.94798839       -0.98443956        1.41551366
 16        -1.68341565       -2.59460979        1.08646957
  6         2.38236576        1.86934843       -0.94929788
  7         3.42670617        2.26016826       -1.27763612
  6         0.68637086        3.24032839        0.86396510
  7         0.96152708        4.25733674        1.35479012
  6         0.17031748       -3.32486427       -0.83470863
  7         0.69156839       -4.18893845       -1.41157016
  6        -0.74312428       -0.54647433        5.31887359
  6        -0.04340914        0.29590936        4.40968985
 16        -0.04907756       -1.28207955        6.75612818
 16         1.65172121        0.77012700        4.60005542
  6         1.12018578       -0.03390492        7.22410835
  6         1.79402412        0.77743660        6.36660565
  6        -2.02879881       -0.75564675        4.86215581
  7        -0.69470504        0.71385166        3.35186287
 16        -2.23361152        0.07929811        3.35691649
  6         1.34230402        0.02545245        8.63060394
  7         1.50756015        0.01985101        9.78119164
  6         2.75158337        1.71971008        6.84529439
  7         3.56075344        2.48048933        7.18790913
  6        -3.05454090       -1.53526047        5.44971351
  7        -3.90501535       -2.17805317        5.91297009

@aug-cc-pVTZ.gbs/N

