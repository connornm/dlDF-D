%Chk=BLIND00831
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.35139641       -0.91175043        0.29383776
  6         0.80265590       -0.48558088        1.00963761
 16        -1.57441056        0.15421376       -0.38179039
 16         1.21366893        1.20998879        1.31080223
  6        -0.57868502        1.55876894       -0.80591002
  6         0.52564453        1.97566320       -0.13189261
  6        -0.40139895       -2.29101111        0.27279204
  7         1.58524057       -1.41544077        1.50033287
 16         0.95443630       -2.91632685        1.15342679
  6        -1.07405829        2.27571980       -1.93359635
  7        -1.52950700        2.84028128       -2.84190787
  6         1.22997302        3.15019745       -0.52980293
  7         1.81708468        4.11613398       -0.80034432
  6        -1.38038966       -3.12632169       -0.31797176
  7        -2.17932892       -3.82702280       -0.78893878
  6         0.88395456       -0.46854949        6.77322588
  6        -0.25222796       -1.13591344        7.31131630
 16         0.89006626        0.49479502        5.30338009
 16        -1.87006213       -1.07594379        6.59465418
  6        -0.75166724        1.16315256        5.34954718
  6        -1.84443406        0.53721011        5.86111415
  6         1.99530999       -0.75137549        7.54147179
  7        -0.07083713       -1.86415559        8.38584911
 16         1.53483959       -1.82397425        8.82300426
  6        -0.84997870        2.44797547        4.74073309
  7        -0.88479643        3.48317588        4.21315289
  6        -3.13319740        1.14545863        5.80575404
  7        -4.20596288        1.59004456        5.75632009
  6         3.32708583       -0.30642337        7.35828021
  7         4.42719699        0.04162430        7.21787301

@aug-cc-pVTZ.gbs/N

