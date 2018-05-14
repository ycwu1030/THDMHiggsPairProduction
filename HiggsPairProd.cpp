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
    cout<<FTriFermion(shat,MT)<<"  "<<FTriFermion(shat,MB)<<endl;
    cout<<FBoxFermion(shat,-500*500,-500*500,MT,Mh,Mh)<<"  "<<FBoxFermion(shat,-500*500,-500*500,MB,Mh,Mh)<<endl;
    cout<<GBoxFermion(shat,-500*500,-500*500,MT,Mh,Mh)<<"  "<<GBoxFermion(shat,-500*500,-500*500,MB,Mh,Mh)<<endl;
    sp.Get_CS_Parton(shat,50,1,1);
    for (double pt = 50; pt < sqrt(shat/4.0-Mh*Mh); pt+=1)
    {
        outfile<<pt<<"  "<<sp.Get_CS_Parton(shat,pt,1,1)<<endl;
    }
    ltexi();
    return 0;
}