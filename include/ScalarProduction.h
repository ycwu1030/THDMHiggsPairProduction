#ifndef ScalarProduction_H
#define ScalarProduction_H

#include "THDMModelParameters.h"
#include "Loopfuncs.h"
#include "cubature.h"
#include "cuba.h" // The Cuba Library
// Using 2HDMC for calculating the model parameters
#include "THDM.h"
#include "DecayTable.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_integration.h>
#include "LHAPDF/LHAPDF.h"

enum CUBAINTEGRATOR
{
    VEGAS = 0,
    SUAVE = 1,
    DIVONNE = 2,
    CUHRE = 3
};

class ScalarProduction
{
public:
    ScalarProduction();
    ~ScalarProduction();

    void Set_Yukawa_Type(int type);
    bool Set_THDM_Params(double mH, double mA, double mHc, double alpha, double beta, double m122); // For Physics basis, not considering lambda6 and lambda7
    //Note that in 2HDMC, the convention is -1<=sba<=1, and cba>=0;
    void Set_PDFSet(std::string pdfname);
    double CS_pp2SS_MISER(double s, int H1, int H2);
    double CS_pp2SS_VEGAS(double s, int H1, int H2);
    double CS_pp2SS_STEPBYSTEP(double s, int H1, int H2);
    double CS_pp2SS_HCUBATURE(double s,int H1, int H2);
    double CS_pp2SS_PCUBATURE(double s,int H1, int H2);
    double CS_pp2SS_CUBAVEGAS(double s,int H1, int H2, CUBAINTEGRATOR Choice);


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
    const double Prefactor = _GF*_GF/(256.0*pow(2*PI,3));
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
    LHAPDF::AlphaS_ODE as_ode;

    double dSigmahatgg2SSdthatGeneral(double shat, double that, int H1, int H2);
    double dSigmahatgg2SSdptGeneral(double shat, double pt, int H1, int H2);
    double PARTONINTEGRAL(double mhh, double S);
    double Sigmagg2SSGeneral(double shat,int H1,int H2);
};

#endif