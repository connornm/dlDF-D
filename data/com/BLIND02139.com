%Chk=BLIND02139
%Mem=2GB
%NProcShared=16
#T M05/Gen test Massage Symmetry=None SCF=(tight,maxcyc=60)

M05 opt

0 1
  6        -0.50151015       -0.22696273       -0.85911838
  6         0.27868614       -1.30336236       -0.35070495
 16        -0.11874094        1.47789971       -0.67084084
 16         1.76229750       -1.09468205        0.59266041
  6         0.66909776        1.45757288        0.91756059
  6         1.41550925        0.43637417        1.41522227
  6        -1.56526940       -0.72633445       -1.58323672
  7        -0.11712290       -2.52061707       -0.63263793
 16        -1.48053479       -2.45770419       -1.58556336
  6         0.48794950        2.66817349        1.64743806
  7         0.33260447        3.68239892        2.19370513
  6         2.04721437        0.54061378        2.68959454
  7         2.60341289        0.59317026        3.70881343
  6        -2.57730153       -0.00950293       -2.26669263
  7        -3.41046664        0.56262116       -2.84080573
  6        -0.50151015        0.22696266        4.58017160
  6         0.27868618        1.30336232        5.08858492
 16        -0.11874098       -1.47789977        4.76844931
 16         1.76229753        1.09468206        6.03195031
  6         0.66909771       -1.45757281        6.35685074
  6         1.41550923       -0.43637408        6.85451232
  6        -1.56526938        0.72633434        3.85605320
  7        -0.11712282        2.52061701        4.80665182
 16        -1.48053471        2.45770408        3.85372640
  6         0.48794942       -2.66817335        7.08672832
  7         0.33260436       -3.68239872        7.63299549
  6         2.04721435       -0.54061358        8.12888459
  7         2.60341287       -0.59316997        9.14810350
  6        -2.57730153        0.00950278        3.17259736
  7        -3.41046665       -0.56262133        2.59848432

@aug-cc-pVTZ.gbs/N

