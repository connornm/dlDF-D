%Chk=BLIND02514
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.14490053        0.80218507       -0.61368711
  6         0.23992706       -0.38265383       -1.30208465
 16        -1.07619457        0.85023617        0.87590986
 16        -0.13506189       -2.02062291       -0.74438498
  6        -0.51261767       -0.63892054        1.65624020
  6        -0.14136891       -1.77621040        1.01095716
  6         0.26256647        1.90883403       -1.33104472
  7         0.88273281       -0.23969678       -2.43532076
 16         1.04877837        1.37980546       -2.78253574
  6        -0.52173001       -0.56281196        3.07934037
  7        -0.56917183       -0.46114242        4.23632162
  6         0.25266531       -2.93845582        1.73750003
  7         0.55246600       -3.91700420        2.28839817
  6         0.07682293        3.27713525       -1.01732581
  7        -0.07473646        4.40494341       -0.78031437
  6         0.31884520       -0.62021484        9.63183511
  6        -1.09779749       -0.54806110        9.51460048
 16         1.45511273       -0.43680390        8.30376389
 16        -1.95298331       -0.22606394        7.99816934
  6         0.57883812        0.72273720        7.28802311
  6        -0.77275298        0.80291621        7.16795280
  6         0.66181692       -0.92699082       10.93320444
  7        -1.80426218       -0.76806870       10.59643426
 16        -0.78655486       -1.11979689       11.86587350
  6         1.44187683        1.57403457        6.53865895
  7         2.18566090        2.22860112        5.93077663
  6        -1.37769113        1.74438607        6.28385651
  7        -1.91233124        2.47512276        5.55508358
  6         1.95456853       -1.10479280       11.48306641
  7         3.00753184       -1.26335439       11.94913764

@aug-cc-pVTZ.gbs/N
