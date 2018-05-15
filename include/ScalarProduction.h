#ifndef ScalarProduction_H
#define ScalarProduction_H

#include "THDMModelParameters.h"
#include "Loopfuncs.h"
// Using 2HDMC for calculating the model parameters
#include "THDM.h"
#include "DecayTable.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include "LHAPDF/LHAPDF.h"


class ScalarProduction
{
public:
    ScalarProduction();
    ~ScalarProduction();

    void Set_Yukawa_Type(int type);
    bool Set_THDM_Params(double mH, double mA, double mHc, double alpha, double beta, double m122); // For Physics basis, not considering lambda6 and lambda7
    //Note that in 2HDMC, the convention is -1<=sba<=1, and cba>=0;
    void Set_PDFSet(std::string pdfname);
    double CS_pp2SS_PT(double sqrts, int H1, int H2);
    double CS_pp2SS_Ctheta(double sqrts, int H1, int H2); // This returns the CS for the process p p > h1 h2; The index here can be used to extend to more general cases.
    // 1,2,3,4 for h,H,A,H+ // All constant will be added in this function, e.g. pre-factor and unit transfer
    // double Get_dSigmapp2SSdMhhdptdxGeneral(double Mhh, double pt, double x, double s, double H1, double H2);
    // double Get_dSigmapp2SSdMhhdCthetadxGeneral(double Mhh, double Ctheta, double x, double s, double H1, double H2);
    double Get_CS_Partial_Ctheta(double *X, size_t dim, void * index);

    double Get_CS_Parton_dt(double shat, double that){return dSigmahatgg2SSdthatGeneral(shat,that,1,1);}
// private:
    THDM _mod; // The THDM model object which will provide some calculations in the 2HDM model
    DecayTable* _dt; // The Decay table of the 2HDM scalars, the width will be used in the cross section calculation.
    int _type;
    double _alpha;
    double _beta;
    double _m122;
    double _mHH;
    double _mA;
    double _mHc;
    double ScalarMasses[4];
    double _Gammah;
    double _GammaH;
    double _GammaA;
    double _GammaHc;
    double Gammas[4];
    const double Prefactor = _GF*_GF*Alfas*Alfas/(256.0*pow(2*PI,3));
    const double GeV2tofb = 0.3894*pow(10,12);
    double _eta; // Symmetric factor accounting for whether the final states are identical particle

    double Get_LambdaAndSymmetryFactor(int H1, int H2, int H3);
    void Calc_Lambda();
    void Calc_Lambdahhh();
    void Calc_LambdahhH();
    void Calc_LambdahHH();
    void Calc_LambdaHHH();
    double Lhhh,LhhH,LhHH,LHHH;

    //double Get_Xi(int Hi);
    void Calc_Xi();
    double xiUp[4],xiDown[4],xiLep[4],xiV[4];//1,2,3,4 for h H A Hc

    void CalcWidth();

    LHAPDF::PDF *pdf;

    double dSigmahatgg2SSdthatGeneral(double shat, double that, int H1, int H2);
    double dSigmahatgg2SSdptGeneral(double shat, double pt, int H1, int H2);// X including shat/mhh and pt of the higgs, this will be used for VEGAS integral.
    double dSigmahatgg2SSdCthetaGeneral(double shat, double Ctheta, int H1, int H2);
    double dSigmapp2SSdMhhdptdxGeneral(double Mhh, double pt, double x, double s, int H1, int H2); //The same as above, but adding the PDF functions using lhapdf
    double dSigmapp2SSdMhhdCthetadxGeneral(double Mhh, double Ctheta, double x, double s, int H1, int H2);
    double dSigmapp2SSdx1dx2dptGeneral(double x1, double x2, double pt, double s, int H1, int H2);
    double dSigmapp2SSdx1dx2dCthetaGeneral(double x1, double x2, double Ctheta, double s, int H1, int H2);
};

#endif