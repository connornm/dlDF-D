%Chk=CH400152
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage SCF=(tight,maxcyc=40)

M05 opt

0 1
  6-Bq        0.00000000       0.00000000       0.00000000
  1-Bq        0.69070410       0.00000000      -0.86539806
  1-Bq       -0.31235745      -1.03861687       0.22292101
  1-Bq        0.51383291       0.42831271       0.88233209
  1-Bq       -0.89217955       0.61030416      -0.23985504
  6         0.00000000        0.00000000        4.43773029
  1        -0.25556813       -0.89645157        5.03526311
  1         0.55250873        0.72064780        5.07128124
  1         0.63254831       -0.29695899        3.57884499
  1        -0.92948891        0.47276276        4.06553181

@aug-cc-pVTZ.gbs/N

