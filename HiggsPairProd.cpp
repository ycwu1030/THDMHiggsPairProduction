#include <iostream>
#include <fstream>
#include "THDMModelParameters.h"
#include "Loopfuncs.h"
#include "ScalarProduction.h"

using namespace std;

int main(int argc, char const *argv[])
{
    ofstream outfile("THDMtest.dat");
    ltini();
    ScalarProduction sp;
    sp.Set_Yukawa_Type(1);
    double shat = 250*250;
    cout<<sp.CS_pp2SS_PT(14000,1,1)<<endl;
    ltexi();
    return 0;
}