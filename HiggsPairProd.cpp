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
    double shat = 1000.0*1000.0;
    for (double pt = 50; pt < sqrt(shat/4.0-Mh*Mh); pt+=1)
    {
        outfile<<pt<<"  "<<sp.Get_CS_Parton(shat,pt,1,1)<<endl;
    }
    ltexi();
    return 0;
}