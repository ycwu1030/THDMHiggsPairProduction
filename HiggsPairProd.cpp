#include <iostream>
#include <fstream>
#include <time.h>
#include "THDMModelParameters.h"
#include "Loopfuncs.h"
#include "ScalarProduction.h"

using namespace std;

int main(int argc, char const *argv[])
{
    clock_t T1, T2, T3;
    time_t start,end;
    ltini();
    ScalarProduction sp;
    sp.Set_Yukawa_Type(1);
    double mphi = 200;
    double cba = 0.1;
    double tb = 1;
    double beta = atan(tb);
    double alpha = beta - acos(cba);
    for (double m12 = 0; m12<= 200; m12+=20)
    {
    cout<<"<<<<<<<<<<  m122 = "<<m12<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
        sp.Set_THDM_Params(mphi,mphi,mphi,alpha,beta,m12*m12);
    cout<<"Cross Section calculation Using Pt of Higgs: "<<endl;
    // T1=clock();
    // cout<<"  STEP BY STEP: "<<sp.CS_pp2SS_STEPBYSTEP(14000*14000,1,1)<<endl;
    // T2=clock();
    // cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    T1=clock(); time(&start);
    cout<<"  GSL MISER: "<<sp.CS_pp2SS_MISER(14000*14000,1,1)<<endl;
    T2=clock(); time(&end);
    cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<"  "<<end-start<<endl;
    // T1=clock();
    // cout<<"  GSL VEGAS: "<<sp.CS_pp2SS_VEGAS(14000*14000,1,1)<<endl;
    // T2=clock();
    // cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    // T1=clock();
    // cout<<"  HCUBATURE: "<<sp.CS_pp2SS_HCUBATURE(14000*14000,1,1)<<endl;
    // T2=clock();
    // cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    // T1=clock();
    // cout<<"  PCUBATURE: "<<sp.CS_pp2SS_PCUBATURE(14000*14000,1,1)<<endl;
    // T2=clock();
    // cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    // T1=clock();
    // cout<<"  CUBA VEGAS: "<<sp.CS_pp2SS_CUBA(14000*14000,1,1,VEGAS)<<endl;
    // T2=clock();
    // cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    // T1=clock();
    // cout<<"  CUBA SUAVE: "<<sp.CS_pp2SS_CUBA(14000*14000,1,1,SUAVE)<<endl;
    // T2=clock();
    // cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    // T1=clock();
    // cout<<"  CUBA DIVONNE: "<<sp.CS_pp2SS_CUBA(14000*14000,1,1,DIVONNE)<<endl;
    // T2=clock();
    // cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    // T1=clock();
    // cout<<"  CUBA CUHRE: "<<sp.CS_pp2SS_CUBA(14000*14000,1,1,CUHRE)<<endl;
    // T2=clock();
    // cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    cout<<"Cross Section calculation Using ctheta: "<<endl;
    // T1=clock();
    // cout<<"  STEP BY STEP: "<<sp.CS_pp2SS_STEPBYSTEPFROMCTHETA(14000*14000,1,1)<<endl;
    // T2=clock();
    // cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    T1=clock(); time(&start);
    cout<<"  GSL MISER: "<<sp.CS_pp2SS_MISERFROMCTHETA(14000*14000,1,1)<<endl;
    T2=clock(); time(&end);
    cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<"  "<<end-start<<endl;
    // T1=clock();
    // cout<<"  GSL VEGAS: "<<sp.CS_pp2SS_VEGASFROMCTHETA(14000*14000,1,1)<<endl;
    // T2=clock();
    // cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    T1=clock(); time(&start);
    cout<<"  HCUBATURE: "<<sp.CS_pp2SS_HCUBATUREFROMCTHETA(14000*14000,1,1)<<endl;
    T2=clock(); time(&end);
    cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<"  "<<end-start<<endl;
    // T1=clock();
    // cout<<"  PCUBATURE: "<<sp.CS_pp2SS_PCUBATUREFROMCTHETA(14000*14000,1,1)<<endl;
    // T2=clock();
    // cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    T1=clock(); time(&start);
    cout<<"  CUBA VEGAS: "<<sp.CS_pp2SS_CUBAFROMCTHETA(14000*14000,1,1,VEGAS)<<endl;
    T2=clock(); time(&end);
    cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<"  "<<end-start<<endl;
    T1=clock(); time(&start);
    cout<<"  CUBA SUAVE: "<<sp.CS_pp2SS_CUBAFROMCTHETA(14000*14000,1,1,SUAVE)<<endl;
    T2=clock(); time(&end);
    cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<"  "<<end-start<<endl;
    // T1=clock();
    // cout<<"  CUBA DIVONNE: "<<sp.CS_pp2SS_CUBAFROMCTHETA(14000*14000,1,1,DIVONNE)<<endl;
    // T2=clock();
    // cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
    T1=clock(); time(&start);
    cout<<"  CUBA CUHRE: "<<sp.CS_pp2SS_CUBAFROMCTHETA(14000*14000,1,1,CUHRE)<<endl;
    T2=clock(); time(&end);
    cout<<"       TIME CONSUMING: "<<(double)(T2-T1)/CLOCKS_PER_SEC<<"  "<<end-start<<endl;
    cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
}
    ltexi();
    return 0;
}
