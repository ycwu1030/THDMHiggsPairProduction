#include <iostream>
#include "ScalarProduction.h"

double lambda(double a,double b,double c)
{
    return a*a+b*b+c*c-2*a*b-2*b*c-2*c*a;
}
ScalarProduction::ScalarProduction()
{
    _type = 2; // Default working in Type-II
    _beta = PI/4.0;
    _alpha = _beta-PI/2.0;
    _m122 = 200.0*200.0;
    _mHH = 200.0;
    _mA = 200.0;
    _mHc = 200.0;
    ScalarMasses[0] = Mh;
    ScalarMasses[1] = _mHH;
    ScalarMasses[2] = _mA;
    ScalarMasses[3] = _mHc;
    _mod.set_param_phys(Mh, _mHH, _mA, _mHc,
                      sin(_beta-_alpha),0,0,
                      _m122, tan(_beta));
    _mod.set_yukawas_type(_type);
    _dt = new DecayTable(_mod);
    CalcWidth();
    Calc_Lambda();
    Calc_Xi();
    pdf = LHAPDF::mkPDF("MSTW2008lo68cl");
}
ScalarProduction::~ScalarProduction()
{ 
    delete _dt;
    delete pdf;
    _dt = NULL;
    pdf = NULL;
}

void ScalarProduction::CalcWidth()
{
    _dt->set_model(_mod);
    _Gammah = _dt->get_gammatot_h(1);
    _GammaH = _dt->get_gammatot_h(2);
    // _GammaA = _dt.get_gammatot_h(3); // For now, didn't using it
    // _GammaHc = _dt.get_gammatot_h(4); // For now, didn't using it
    Gammas[0] = _Gammah;
    Gammas[1] = _GammaH;
    // Gammas[2] = _GammaA;
    // Gammas[3] = _GammaHc;
}

void ScalarProduction::Set_Yukawa_Type(int type)
{
    _type = type;
    _mod.set_yukawas_type(_type);
}

bool ScalarProduction::Set_THDM_Params(double mH, double mA, double mHc, double alpha, double beta, double m122)
{
    _beta = beta;
    _alpha = alpha;
    _m122 = m122;
    _mHH = mH;
    _mA = mA;
    _mHc = mHc;
    ScalarMasses[1] = _mHH;
    ScalarMasses[2] = _mA;
    ScalarMasses[3] = _mHc;
    bool pset = _mod.set_param_phys(Mh, _mHH, _mA, _mHc,
                      sin(_beta-_alpha),0,0,
                      _m122, tan(_beta));
    if (!pset)
    {
        return false;
    }
    _mod.set_yukawas_type(_type);
    CalcWidth();
    Calc_Lambda();
    Calc_Xi();
    return true;
}

void ScalarProduction::Set_PDFSet(std::string pdfname)
{
    delete pdf;
    pdf = LHAPDF::mkPDF(pdfname);
}

double ScalarProduction::Get_LambdaAndSymmetryFactor(int H1, int H2, int H3)
{
    // 1,2,3,4 for h H A Hc
    // Implemented:
    // hhh vertex = 1*1*1 = 1;
    // hhH vertex = 1*1*2 = 2;
    // hHH vertex = 1*2*2 = 4;
    // HHH vertex = 2*2*2 = 8; 
    int vertex = H1*H2*H3;
    switch(vertex){
        case 1: return 3.0*Lhhh;
        case 2: return 2.0*LhhH;
        case 4: return 2.0*LhHH;
        case 8: return 3.0*LHHH;
        default: return 3.0*Lhhh;
    }
}

void ScalarProduction::Calc_Lambda()
{
    Calc_Lambdahhh();
    Calc_LambdahhH();
    Calc_LambdahHH();
    Calc_LambdaHHH();
}

void ScalarProduction::Calc_Lambdahhh()
{
    double sb = sin(_beta);
    double cb = cos(_beta);
    double sab = sin(_alpha-_beta);
    double s3a3b = sin(3*_alpha-3*_beta);
    double s3a1b = sin(3*_alpha+_beta);
    double s1a3b = sin(_alpha+3*_beta);
    double c3a1b = cos(3*_alpha-_beta);
    double c1a3b = cos(_alpha-3*_beta);
    double cab = cos(_alpha+_beta);
    double pref = - 1.0/(32.0*vev*sb*sb*cb*cb);
    Lhhh = pref*(Mh*Mh*(3.0*sab + s3a3b - s3a1b - 3.0*s1a3b)+4*_m122*(c3a1b + c1a3b + 2.0*cab));
}
void ScalarProduction::Calc_LambdahhH()
{
    double cab = cos(_alpha - _beta);
    double sb = sin(_beta);
    double cb = cos(_beta);
    double sa = sin(_alpha);
    double ca = cos(_alpha);
    double s2a = sin(2*_alpha);
    double s2b = sin(2*_beta);
    LhhH = cab/(2.0*vev*sb*cb)*((2*Mh*Mh+_mHH*_mHH)*sa*ca + _m122*(1-3*s2a/s2b));
}
void ScalarProduction::Calc_LambdahHH()
{
    double sba = cos(_beta - _alpha);
    double sb = sin(_beta);
    double cb = cos(_beta);
    double sa = sin(_alpha);
    double ca = cos(_alpha);
    double s2a = sin(2*_alpha);
    double s2b = sin(2*_beta);
    LhHH = sba/(2.0*vev*sb*cb)*(-(Mh*Mh+2*_mHH*_mHH)*sa*ca + _m122*(1+3*s2a/s2b));
}
void ScalarProduction::Calc_LambdaHHH()
{
    double sb = sin(_beta);
    double cb = cos(_beta);
    double c3a3b = cos(3*_alpha - 3*_beta);
    double c3a1b = cos(3*_alpha + _beta);
    double cab = cos(_alpha - _beta);
    double c1a3b = cos(_alpha + 3*_beta);
    double s1a3b = sin(_alpha - 3*_beta);
    double s3a1b = sin(3*_alpha - _beta);
    double sab = sin(_alpha + _beta);
    double pref = -1.0/(32.0*vev*sb*sb*cb*cb);
    LHHH = pref*(_mHH*_mHH*(c3a3b - c3a1b - 3*cab + 3*c1a3b)+4*_m122*(s1a3b-s3a1b+2*sab));
}

void ScalarProduction::Calc_Xi()
{
// Type I and Type II implemented
    switch(_type){
        case 1: { 
            xiUp[0] = cos(_alpha)/sin(_beta); 
            xiDown[0] = cos(_alpha)/sin(_beta); 
            xiLep[0] = cos(_alpha)/sin(_beta); 
            xiV[0] = sin(_beta-_alpha);
            xiUp[1] = sin(_alpha)/sin(_beta);
            xiDown[1] = sin(_alpha)/sin(_beta);
            xiLep[1] = sin(_alpha)/sin(_beta);
            xiV[1] = cos(_beta-_alpha);
            return;
        }
        case 2: {
            xiUp[0] = cos(_alpha)/sin(_beta);
            xiDown[0] = -sin(_alpha)/cos(_beta);
            xiLep[0] = -sin(_alpha)/cos(_beta);
            xiV[0] = sin(_beta-_alpha);
            xiUp[1] = sin(_alpha)/sin(_beta);
            xiDown[1] = cos(_alpha)/cos(_beta);
            xiLep[1] = cos(_alpha)/cos(_beta);
            xiV[1] = cos(_beta-_alpha);
            return;
        }
        default: { // Default return Type-II
            xiUp[0] = cos(_alpha)/sin(_beta);
            xiDown[0] = -sin(_alpha)/cos(_beta);
            xiLep[0] = -sin(_alpha)/cos(_beta);
            xiV[0] = sin(_beta-_alpha);
            xiUp[1] = sin(_alpha)/sin(_beta);
            xiDown[1] = cos(_alpha)/cos(_beta);
            xiLep[1] = cos(_alpha)/cos(_beta);
            xiV[1] = cos(_beta-_alpha);
            return;
        }
    }
}

double ScalarProduction::dSigmahatgg2SSdthatGeneral(double shat, double that, int H1, int H2)
{
    complex<double> img(0.0,1.0);
    double mc = ScalarMasses[H1-1];
    double md = ScalarMasses[H2-1];
    double lam = lambda(shat,mc*mc,md*md);
    double uhat = (pow(mc,2)+pow(md,2))-shat-that;
    double c1 = Get_LambdaAndSymmetryFactor(H1,H2,1);// The symmetry factor is included here
    double c2 = Get_LambdaAndSymmetryFactor(H1,H2,2);
#ifdef DEBUG
    std::cout<<Jacobi<<std::endl;
    std::cout<<c1<<"  "<<c2<<std::endl;
#endif
    ComplexType CTriTop = c1*vev*xiUp[0]/(shat - Mh * Mh + img*Mh*_Gammah) + c2*vev*xiUp[1]/(shat - _mHH*_mHH + img *_mHH*_GammaH);
    ComplexType CTriBot = c1*vev*xiDown[0]/(shat - Mh * Mh + img*Mh*_Gammah) + c2*vev*xiDown[1]/(shat - _mHH*_mHH + img *_mHH*_GammaH);
    ComplexType CBoxTop = xiUp[H1-1]*xiUp[H2-1];
    ComplexType CBoxBot = xiDown[H1-1]*xiDown[H2-1];
    ComplexType P1 = CTriTop*FTriFermion(shat,MT) + CTriBot*FTriFermion(shat,MB) + CBoxTop*FBoxFermion(shat,that,uhat,MT,mc,md) + CBoxBot*FBoxFermion(shat,that,uhat,MB,mc,md);
    ComplexType P2 = CBoxTop*GBoxFermion(shat,that,uhat,MT,mc,md) + CBoxTop*GBoxFermion(shat,that,uhat,MB,mc,md);
    return (norm(P1) + norm(P2));
}

double ScalarProduction::dSigmahatgg2SSdptGeneral(double shat, double pt, int H1, int H2)
{
    double mc = ScalarMasses[H1-1];
    double md = ScalarMasses[H2-1];
    double lam = lambda(shat,mc*mc,md*md);
    if (lam < 4*shat*pt*pt) // Only consider the on-shell production
    {
        return 0.0;
    }
    double that = (pow(mc,2)+pow(md,2))/2.0 - shat/2 + 1.0/2.0*sqrt(lam-4*shat*pt*pt);
    double Jacobi = 2.0*pt*shat/(sqrt(lam-4*shat*pt*pt));
    return Jacobi*dSigmahatgg2SSdthatGeneral(shat,that,H1,H2);
}

double ScalarProduction::dSigmahatgg2SSdCthetaGeneral(double shat, double Ctheta, int H1, int H2)
{
// X[0] is the shat i.e. Mhh2;
// X[1] is the cos\theta polar angle of the Higgs ;
    double mc = ScalarMasses[H1-1];
    double md = ScalarMasses[H2-1];
    double lam = lambda(shat,mc*mc,md*md);
    double that = (pow(mc,2)+pow(md,2))/2.0 - shat/2 + 1.0/2.0*sqrt(lam)*Ctheta;
    double Jacobi = 1.0/2.0*sqrt(lam);
    return Jacobi*dSigmahatgg2SSdthatGeneral(shat,that,H1,H2);
}

double ScalarProduction::dSigmapp2SSdMhhdptdxGeneral(double Mhh, double pt, double x, double s, int H1, int H2)
{
    double shat = Mhh*Mhh;
    double tau = shat/s;
    if (x<tau || x>1)
    {
        return 0;
    }
    double pdf1 = pdf->xfxQ2(21,x,shat); // this is x*f1
    double pdf2 = pdf->xfxQ2(21,tau/x,shat); // this is tau/x*f2
    double partoncs = dSigmahatgg2SSdptGeneral(shat,pt,H1,H2);
    return 1.0/x*pdf1*pdf2*2.0/Mhh*partoncs;
}

double ScalarProduction::dSigmapp2SSdMhhdCthetadxGeneral(double Mhh, double Ctheta, double x, double s, int H1, int H2)
{
    double shat = Mhh*Mhh;
    double tau = shat/s;
    if (x<tau || x>1)
    {
        return 0;
    }
    double pdf1 = pdf->xfxQ2(21,x,shat); // this is x*f1
    double pdf2 = pdf->xfxQ2(21,tau/x,shat); // this is tau/x*f2
    double partoncs = dSigmahatgg2SSdCthetaGeneral(shat,Ctheta,H1,H2);
    return 1.0/x*pdf1*pdf2*2.0/Mhh*partoncs;
}

typedef struct{
    int H1;
    int H2;
    double s;
    ScalarProduction *sp;
}VOIDVEGAS;

double partialCS_for_VEGAS_PT(double *X, size_t dim, void *modparameter)
{
    VOIDVEGAS *fp = (VOIDVEGAS*)modparameter;
    return (fp->sp)->dSigmapp2SSdMhhdptdxGeneral(X[0],X[1],X[2],fp->s,fp->H1,fp->H2);
}
double partialCS_for_VEGAS_Ctheta(double *X, size_t dim, void *modparameter)
{
    VOIDVEGAS *fp = (VOIDVEGAS*)modparameter;
    return (fp->sp)->dSigmapp2SSdMhhdCthetadxGeneral(X[0],X[1],X[2],fp->s,fp->H1,fp->H2);
}

double ScalarProduction::CS_pp2SS_PT(double sqrts, int H1, int H2)
{
    _eta = H1==H2?1.0/2.0:1.0;
    double res, err;

    double Q2L[3] = {ScalarMasses[H1-1]+ScalarMasses[H2-1]+0.1,0,0};
    double Q2U[3] = {sqrts-0.1,sqrts/2.0-0.1,1};

    const gsl_rng_type *T;
    gsl_rng *r;

    VOIDVEGAS fp = {H1,H2,sqrts*sqrts,this};
    gsl_monte_function G = {&partialCS_for_VEGAS_PT, 3, &fp};

    size_t calls = 10000;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(3);

    gsl_monte_vegas_integrate(&G,Q2L,Q2U,3,calls,r,s,&res,&err);

    do{
        gsl_monte_vegas_integrate(&G,Q2L,Q2U,3,calls,r,s,&res,&err);
    }while(fabs(gsl_monte_vegas_chisq(s)-1.0)>0.5);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);
    s=NULL;
    r=NULL;

    return res*_eta*GeV2tofb*Prefactor;
}

double ScalarProduction::CS_pp2SS_Ctheta(double sqrts, int H1, int H2)
{
    _eta = H1==H2?1.0/2.0:1.0;
    double res, err;

    double Q2L[3] = {ScalarMasses[H1-1]+ScalarMasses[H2-1]+0.1,-1,0};
    double Q2U[3] = {sqrts-0.1,1,1};

    const gsl_rng_type *T;
    gsl_rng *r;

    VOIDVEGAS fp = {H1,H2,sqrts*sqrts,this};
    gsl_monte_function G = {&partialCS_for_VEGAS_Ctheta, 3, &fp};

    size_t calls = 1000;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(3);

    gsl_monte_vegas_integrate(&G,Q2L,Q2U,3,calls,r,s,&res,&err);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);
    s=NULL;
    r=NULL;

    return res*_eta*GeV2tofb*Prefactor;
}
