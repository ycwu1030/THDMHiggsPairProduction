#ifndef THDMParameter_H
#define THDMParameter_H

#include <cmath>

#define DEBUG 1

#define PI 3.1415926535897932384626433832795029
#define PI2 (PI*PI)
#define PiHalf (PI/2)

#define MZ 91.1876
#define GAZ 2.4952
#define MZ2 (MZ*MZ)
#define MW 80.385
#define GAW 2.085
#define MW2 (MW*MW)
#define cw (MW/MZ)
#define w (acos(cw))
#define sw (sin(w))
#define cw2 (cw*cw)
#define sw2 (sw*sw)

#define GF 1.16637e-5
#define Alfa0 (1/137.035999074)
#define AlfaGF (sqrt2/pi*GF*MW2*SW2)
#define AlfaMZ (1/127.944)
#define Alfa Alfa0
#define Alfa2 (Alfa*Alfa)
#define DeltaAlfa5Had .027547
#define AlfasMZ .1184
#define Alfas 0.0
#define Alfashgg 0.1184

#define v 246.2
#define v2 (v*v)


#define Mh 125.09
#define Mh02 (Mh*Mh)
#define Mh2 (Mh*Mh)

#define ME .5109989280e-3
#define ME2 (ME*ME)
#define MM 105.6583715e-3
#define MM2 (MM*MM)
#define ML 1776.82e-3
#define ML2 (ML*ML)


#define MU 7.356e-2
#define MU2 (MU*MU)
#define MC 1.275
#define MC2 (MC*MC)
#define MT 173.21
#define MT2 (MT*MT)

#define MD MU
#define MD2 (MD*MD)
#define MS 95e-3
#define MS2 (MS*MS)
#define MB1S 4.66
#define MBatMB 4.18
#define MB 4.66
#define MB2 (MB*MB)

#define Sin(i) sin(i)
#define Csc(i) (1/sin(i))
#define Cos(i) cos(i)
#define Sec(i) (1/cos(i))
#define Tan(i) tan(i)
#define Cot(i) (1/tan(i))

#endif