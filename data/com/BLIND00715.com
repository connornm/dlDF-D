%Chk=BLIND00715
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.73643256        0.65725481        0.25844171
  6         0.06626133        1.09706513       -0.83154731
 16        -1.09433550       -1.02056487        0.63968081
 16         0.86977793        0.00959019       -1.97450298
  6         0.40654026       -1.79934764        0.10583494
  6         1.18251178       -1.38918228       -0.93208444
  6        -1.24927164        1.75009558        0.92776848
  7         0.18394769        2.39013356       -1.01038792
 16        -0.71529941        3.19630627        0.13519364
  6         0.72758403       -2.96821435        0.85534889
  7         0.93464182       -3.92713456        1.47888087
  6         2.34922527       -2.11701609       -1.31021076
  7         3.29094316       -2.69463145       -1.67139498
  6        -2.10222373        1.76803553        2.05795401
  7        -2.80971200        1.80152544        2.97960182
  6         0.16886817       -0.82038298        6.73273669
  6        -1.18840638       -0.39189101        6.72753161
 16         1.29364233       -0.63622769        5.39501977
 16        -1.96288854        0.43481564        5.36689873
  6         0.71250124        0.88466308        4.69267486
  6        -0.57979753        1.30596196        4.68194135
  6         0.44712668       -1.46616695        7.92052391
  7        -1.90958140       -0.66075306        7.78841795
 16        -0.98888978       -1.50221399        8.89077331
  6         1.74584129        1.64536884        4.07249364
  7         2.61759006        2.21795623        3.55927752
  6        -0.94934809        2.52808377        4.04649459
  7        -1.29997849        3.50101690        3.51603669
  6         1.66541241       -2.05649131        8.33581355
  7         2.65464470       -2.55417746        8.68911781

@aug-cc-pVTZ.gbs/N
