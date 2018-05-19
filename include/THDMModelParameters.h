#ifndef THDMParameter_H
#define THDMParameter_H

#include <cmath>

// #define DEBUG 1

#define PI 3.1415926535897932384626433832795029
#define PI2 (PI*PI)
#define PiHalf (PI/2)

#define _MZ 91.1876
#define GAZ 2.4952
#define MZ2 (_MZ*_MZ)
#define _MW 80.385
#define GAW 2.085
#define MW2 (_MW*_MW)
#define cw (_MW/_MZ)
#define thw (acos(cw))
#define sw (sin(thw))
#define cw2 (cw*cw)
#define sw2 (sw*sw)

#define _GF 1.16637e-5
#define Alfa0 (1/137.035999074)
#define AlfaGF (sqrt2/pi*_GF*MW2*SW2)
#define AlfaMZ (1/127.944)
#define Alfa Alfa0
#define Alfa2 (Alfa*Alfa)
#define DeltaAlfa5Had .027547
#define AlfasMZ .1184
#define Alfas 0.1184
#define Alfashgg 0.1184

#define vev 246.2
#define vev2 (vev*vev)


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

#define GSLCALLS 50000

//--- CUBA VARIABLES
#define NDIM 3
#define NCOMP 1
#define NVEC 1
#define EPSREL 1e-2
#define EPSABS 1e-1
#define VERBOSE 0
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 5000000

#define NSTART 2000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 10.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

#endif