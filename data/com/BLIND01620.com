%Chk=BLIND01620
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.11368984        0.05921181       -1.01226355
  6        -0.30474184       -1.25127427       -0.49077336
 16        -0.46598825        1.55162523       -0.15371557
 16        -0.90564188       -1.59289373        1.13919708
  6        -0.10322272        1.05905738        1.51058614
  6        -0.27992310       -0.18843732        2.02066260
  6         0.32241542       -0.02048398       -2.31946055
  7        -0.04809949       -2.26030499       -1.28692142
 16         0.42602349       -1.68926786       -2.77690158
  6         0.36762886        2.13219255        2.32161594
  7         0.72479036        3.04334134        2.94886117
  6        -0.00051039       -0.46881496        3.39079356
  7         0.18613572       -0.73538176        4.50660336
  6         0.62272580        1.03698743       -3.21206143
  7         0.86469796        1.89238107       -3.96095870
  6        -0.11792625        0.05928869        9.38822573
  6        -0.30784769       -1.25103094        9.91054647
 16        -0.46535682        1.55200601       10.24822786
 16        -0.90212172       -1.59212835       11.54305352
  6        -0.09595319        1.05916704       11.91098845
  6        -0.27151959       -0.18817325       12.42183354
  6         0.31257424       -0.02078792        8.07919545
  7        -0.05540994       -2.26028589        9.11333905
 16         0.41286993       -1.68966414        7.62135324
  6         0.37921203        2.13192741       12.71999577
  7         0.73977661        3.04279206       13.34570480
  6         0.01346114       -0.46875829       13.79077472
  7         0.20461064       -0.73546028       14.90578954
  6         0.60997326        1.03642100        7.18530979
  7         0.84947664        1.89160264        6.43537752

@aug-cc-pVTZ.gbs/N

