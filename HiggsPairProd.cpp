#include <iostream>
#include <fstream>
#include "THDMModelParameters.h"
#include "Loopfuncs.h"

using namespace std;

int main(int argc, char const *argv[])
{
    ofstream outfile("THDMtest.dat");
    ltini();
    ComplexType FTriTop;
    double shat;
    for (int i = 0; i < 10; ++i)
    {
        shat = pow(500.0 + i*10.0,2);
        FTriTop = FTriFermion(shat,MT);
        outfile<<shat<<"  "<<FTriTop<<endl;
    }
    ltexi();
    return 0;
}