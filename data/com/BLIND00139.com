%Chk=BLIND00139
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.76691375       -0.24683687        0.62611802
  6        -0.66338180        1.14961531        0.37110996
 16        -0.31019102       -1.51778005       -0.49843719
 16        -0.01978590        1.83984567       -1.12687358
  6         1.04463174       -0.73671738       -1.33421136
  6         1.15624302        0.59491900       -1.58307182
  6        -1.34677185       -0.45120173        1.86195895
  7        -1.10562726        1.96584921        1.29633491
 16        -1.72637190        1.08659134        2.56617686
  6         2.03158962       -1.66080140       -1.78484704
  7         2.79679066       -2.45449083       -2.15324558
  6         2.26836503        1.11672302       -2.30753516
  7         3.14102440        1.58179015       -2.91835782
  6        -1.64531027       -1.67921235        2.50073369
  7        -1.90578701       -2.67770576        3.03576518
  6        -0.93806447       -0.40065068        6.92990580
  6        -0.41375023       -0.09107558        8.21640668
 16        -0.55622405        0.46156950        5.44690517
 16         0.74696282        1.20890233        8.52848196
  6         1.12168152        0.93392442        5.77192177
  6         1.63511086        1.23095656        6.99495336
  6        -1.84632523       -1.43468659        7.03542990
  7        -0.84944303       -0.79731987        9.23092667
 16        -1.97982560       -1.89810685        8.70043617
  6         1.91568527        1.02549230        4.59197981
  7         2.51606479        1.10946115        3.60017001
  6         2.99122087        1.64738529        7.14156403
  7         4.08058281        2.01845764        7.30441980
  6        -2.60633936       -2.04546463        6.00854279
  7        -3.24481303       -2.55208173        5.17981702

@aug-cc-pVTZ.gbs/N

