%Chk=CH400261
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6         0.00000000        0.00000000        0.00000000
  1         1.04391278        0.00000000        0.36909627
  1         0.00001637       -0.00000008       -1.10724250
  1        -0.52196460       -0.90405967        0.36907318
  1        -0.52196455        0.90405976        0.36907305
  0        -0.69077240        0.00000000       -0.24423642
  0        -0.00001083        0.00000006        0.73267860
  0         0.34539163        0.59822955       -0.24422114
  0         0.34539160       -0.59822960       -0.24422105
  6         0.00000000        0.00000000        3.64494639
  1        -0.00000001        0.00000054        4.75218889
  1        -1.04391824       -0.00000018        3.27586555
  1         0.52195912        0.90405954        3.27586512
  1         0.52195912       -0.90405990        3.27586600
  0         0.00000000       -0.00000036        2.91226779
  0         0.69077601        0.00000012        3.88917260
  0        -0.34538801       -0.59822946        3.88917288
  0        -0.34538801        0.59822969        3.88917230

@aug-cc-pVTZ.gbs/N

