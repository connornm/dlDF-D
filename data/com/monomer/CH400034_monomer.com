%Chk=CH400034
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.10723899        0.00000000       -0.00278651
  1        -0.37107326       -0.67987265       -0.79124191
  1        -0.36660111       -0.34610584        0.98579934
  1        -0.36956462        1.02597849       -0.19177092
  6-Bq        0.00000000       0.00000000       3.22043238
  1-Bq        0.11896844      -1.03376662       2.84206876
  1-Bq        0.22047508       0.71833799       2.40718751
  1-Bq       -1.03985802       0.14728774       3.57112822
  1-Bq        0.70041451       0.16814089       4.06134505

@aug-cc-pVTZ.gbs/N

