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
}
ScalarProduction::~ScalarProduction()
{ 
    delete _dt;
    _dt = NULL;
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

double ScalarProduction::dSigmahatgg2SSdptGeneral(double *X, size_t dim, void *modparameter)
{
// X[0] is the shat i.e. Mhh;
// X[1] is the pt of the Higgs;
    complex<double> img(0.0,1.0);
    int *finalindex = (int*)modparameter;
    double mc = ScalarMasses[finalindex[0]-1];
    double md = ScalarMasses[finalindex[1]-1];
    double shat = X[0];
    double lam = lambda(X[0],mc*mc,md*md);
    if (lam < 4*X[0]*X[1]*X[1]) // Only consider the on-shell production
    {
        return 0;
    }
    double that = (pow(mc,2)+pow(md,2))/2.0 - X[0]/2 + 1.0/2.0*sqrt(lam-4*X[0]*X[1]*X[1]);
    double uhat = (pow(mc,2)+pow(md,2))/2.0 - X[0]/2 - 1.0/2.0*sqrt(lam-4*X[0]*X[1]*X[1]);
    double Jacobi = 2.0*X[1]*shat/(sqrt(lam-4*X[0]*X[1]*X[1]));
    double c1 = Get_LambdaAndSymmetryFactor(finalindex[0],finalindex[1],1);// The symmetry factor is included here
    double c2 = Get_LambdaAndSymmetryFactor(finalindex[0],finalindex[1],2);
    ComplexType CTriTop = c1*vev*xiUp[0]/(X[0] - Mh * Mh + img*Mh*_Gammah) + c2*vev*xiUp[1]/(X[0] - _mHH*_mHH + img *_mHH*_GammaH);
    ComplexType CTriBot = c1*vev*xiDown[0]/(X[0] - Mh * Mh + img*Mh*_Gammah) + c2*vev*xiDown[1]/(X[0] - _mHH*_mHH + img *_mHH*_GammaH);
    ComplexType CBoxTop = xiUp[finalindex[0]-1]*xiUp[finalindex[1]-1];
    ComplexType CBoxBot = xiDown[finalindex[0]-1]*xiDown[finalindex[1]-1];
    ComplexType P1 = CTriTop*FTriFermion(shat,MT) + CTriBot*FTriFermion(shat,MB) + CBoxTop*FBoxFermion(shat,that,uhat,MT,mc,md) + CBoxBot*FBoxFermion(shat,that,uhat,MB,mc,md);
    ComplexType P2 = CBoxTop*GBoxFermion(shat,that,uhat,MT,mc,md) + CBoxTop*GBoxFermion(shat,that,uhat,MB,mc,md);
    return Jacobi*(norm(P1) + norm(P2));
}

double ScalarProduction::Get_CS_Parton(double shat, double pt,int H1, int H2)
{
    double X[2] = {shat, pt};
    int index[2] = {H1,H2};
    return dSigmahatgg2SSdptGeneral(X,2,(void*)index);
}