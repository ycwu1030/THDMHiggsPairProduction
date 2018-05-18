#include <iostream>
#include <fstream>
#include <time.h>
#include "THDMModelParameters.h"
#include "Loopfuncs.h"
#include "ScalarProduction.h"

using namespace std;

clock_t T1, T2, T3;
int main(int argc, char const *argv[])
{
    ltini();
    ScalarProduction sp;
    sp.Set_Yukawa_Type(1);
    // sp.Set_PDFSet("xxxxx"); // Setting the PDF set Default is MSTW2008lo68cl
    double cba = 0.1;
    double tb = 1.0;
    double beta = atan(tb);
    double alpha = beta-acos(cba);
    double m12 = 500;   
    //Following calculate the cross section:
    // First argument is the center energy square 
    // The second and third are two indexes of two final scalars 1 for h, 2 for H
    // for (m12 = 0; m12 < 100; m12+=10)
    {
        T1 = clock();
        sp.Set_THDM_Params(300, 300, 300, alpha, beta,0); // Set Parameters mH mA mHc, alpha beta m12^2
        cout<<sp._Gammah<<endl;
        cout<<sp._GammaH<<endl;
        // cout<<sp.CS_pp2SS_STEPBYSTEP(14000*14000,1,1)<<endl;
        cout<<sp.CS_pp2SS_MISER(14000*14000, 1,1)<<endl;
        T2 = clock();
        std::cout<<"Time:  "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    }
    
    ltexi();
    return 0;
}