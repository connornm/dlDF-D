%Chk=CH400177
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         0.78350487        0.00000000       -0.78237209
  1        -0.47588635       -0.99871117        0.04576174
  1         0.45733135        0.23619029        0.98033063
  1        -0.76494986        0.76252088       -0.24372029
  0        -0.51845666        0.00000000        0.51770709
  0         0.31490098        0.66086183       -0.03028122
  0        -0.30262286       -0.15629058       -0.64869916
  0         0.50617854       -0.50457125        0.16127329
  6         0.00000000        0.00000000        3.30971923
  1         0.71106565        0.08132107        2.46487695
  1         0.55950579        0.00847931        4.26516014
  1        -0.70165954        0.85619106        3.28529143
  1        -0.56891189       -0.94599144        3.22354840
  0        -0.47052257       -0.05381135        3.86876376
  0        -0.37023319       -0.00561088        2.67748997
  0         0.46429841       -0.56655418        3.32588346
  0         0.37645735        0.62597641        3.36673973

@aug-cc-pVTZ.gbs/N

