%Chk=BLIND00807
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.90084308       -0.04249812        0.47726782
  6        -0.56349996        1.25771101        0.00664855
 16        -0.40647540       -1.55145182       -0.27584306
 16         0.44791981        1.56780347       -1.41302080
  6         1.17803782       -1.08540258       -0.92100818
  6         1.51335570        0.15210507       -1.37242905
  6        -1.73394026        0.06215626        1.57284975
  7        -1.06549616        2.28173271        0.65258925
 16        -2.03522942        1.73911882        1.89198606
  6         2.10162209       -2.16952922       -0.97345992
  7         2.81109188       -3.08920619       -1.01862522
  6         2.80559996        0.41055969       -1.91763353
  7         3.83534497        0.65902786       -2.39597954
  6        -2.31349072       -0.97807523        2.33914482
  7        -2.80502045       -1.81906798        2.97335743
  6        -0.68622598        0.22012074        6.74847581
  6        -1.20197198        0.30461298        5.42458239
 16         0.59243909       -0.86723998        7.26922425
 16        -0.60886421       -0.66187993        4.06491222
  6         1.57072223       -0.94461908        5.79237062
  6         1.09119368       -0.86420701        4.52302815
  6        -1.39121887        1.07037175        7.57634143
  7        -2.20012943        1.12850727        5.21784605
 16        -2.62539803        1.86100487        6.65085364
  6         2.95704141       -1.15179954        6.04982945
  7         4.07366323       -1.34474915        6.30890318
  6         1.96110642       -0.98580162        3.39945861
  7         2.62792100       -1.10846163        2.45540407
  6        -1.21785683        1.28972066        8.96450572
  7        -1.09579350        1.47533562       10.10544548

@aug-cc-pVTZ.gbs/N

