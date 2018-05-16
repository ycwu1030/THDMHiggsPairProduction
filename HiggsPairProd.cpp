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
    ofstream outfile("THDMtest.dat");
    ltini();
    ScalarProduction sp;
    sp.Set_Yukawa_Type(1);
    double shat = 250*250;
    // for (int mhh = 300; mhh < 2000; mhh+=50)
    {
        // for (int i = 0; i < sqrt(mhh*mhh/4.0-Mh*Mh); i+=5)
        // {
            T1 = clock();
            cout<<"  "<<"  "<<sp.CS_pp2SS_VEGAS(14000*14000, 1,1)<<endl;
            T2 = clock();
            std::cout<<(double)(T2-T1)/CLOCKS_PER_SEC<<endl;
        // }
    }
    ltexi();
    return 0;
}