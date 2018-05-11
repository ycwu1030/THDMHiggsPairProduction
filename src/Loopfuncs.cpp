#include "Loopfuncs.h"

// From the Appendix A of hep-ph/9603205
ComplexType FTriFermion(double shat, double mf)
{
    double mf2 = mf*mf;
    double S = shat/mf2;
    return 2.0/S*(2.0+(4.0-S)*mf2*C0(0,0,shat,mf2,mf2,mf2));
}

ComplexType FBoxFermion(double shat, double that, double uhat, double mf, double mhc, double mhd)
{
    double mf2 = mf*mf;
    double mhc2 = mhc*mhc;
    double mhd2 = mhd*mhd;
    double rhoc = mhc2/mf2;
    double rhod = mhd2/mf2;
    double S = shat/mf2;
    double T = that/mf2;
    double U = uhat/mf2;
    double T1 = T - rhoc;
    double U1 = U - rhoc;
    double T2 = T - rhod;
    double U2 = U - rhod;
    return 2/(S*S)*(4.0*S+8.0*S*mf2*C0(0,0,shat,mf2,mf2,mf2)-2.0*S*(S+rhoc+rhod-8)*mf2*mf2
        *(D0(0,0,mhc2,mhd2,shat,uhat,mf2,mf2,mf2,mf2)+D0(0,0,mhc2,mhd2,shat,that,mf2,mf2,mf2,mf2)+D0(0,mhc2,0,mhd2,that,uhat,mf2,mf2,mf2,mf2))
        +(rhoc+rhod-8.0)*mf2*(T1*C0(0,mhc2,that,mf2,mf2,mf2)+U1*C0(0,mhc2,uhat,mf2,mf2,mf2)
            +U2*C0(0,mhd2,uhat,mf2,mf2,mf2)+T2*C0(0,mhd2,that,mf2,mf2,mf2)-
            (T*U-rhoc*rhod)*mf2*D0(0,mhc2,0,mhd2,that,uhat,mf2,mf2,mf2,mf2)));
}

ComplexType GBoxFermion(double shat, double that, double uhat, double mf, double mhc, double mhd)
{
    double mf2 = mf*mf;
    double mhc2 = mhc*mhc;
    double mhd2 = mhd*mhd;
    double rhoc = mhc2/mf2;
    double rhod = mhd2/mf2;
    double S = shat/mf2;
    double T = that/mf2;
    double U = uhat/mf2;
    double T1 = T - rhoc;
    double U1 = U - rhoc;
    double T2 = T - rhod;
    double U2 = U - rhod;
    return 1.0/(S*(T*U-rhoc*rhod))*((T*T+rhoc*rhod-8.0*T)*mf2*(S*C0(0,0,shat,mf2,mf2,mf2)+T1*C0(0,mhc2,that,mf2,mf2,mf2)+T2*C0(0,mhd2,that,mf2,mf2,mf2)-S*T*mf2*D0(0,0,mhc2,mhd2,shat,that,mf2,mf2,mf2,mf2))
        +(U*U+rhoc*rhod-8.0*U)*mf2*(S*C0(0,0,shat,mf2,mf2,mf2)+U1*C0(0,mhc2,uhat,mf2,mf2,mf2)+U2*C0(0,mhd2,uhat,mf2,mf2,mf2)-S*U*mf2*D0(0,0,mhc2,mhd2,shat,uhat,mf2,mf2,mf2,mf2))
        -(T*T+U*U-2.0*rhoc*rhod)*(T+U-8.0)*mf2*C0(mhc2,mhd2,shat,mf2,mf2,mf2)
        -2.0*(T+U-8.0)*(T*U-rhoc*rhod)*mf2*mf2*(D0(0,0,mhc2,mhd2,shat,uhat,mf2,mf2,mf2,mf2)+D0(0,0,mhc2,mhd2,shat,that,mf2,mf2,mf2,mf2)+D0(0,mhc2,0,mhd2,that,uhat,mf2,mf2,mf2,mf2)));
}
