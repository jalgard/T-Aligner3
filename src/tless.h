// by J@:)
#ifndef TLESS_H
#define TLESS_H

#include <vector>
#include <string>
#include <iostream>

using namespace std;



struct TlessDNA
{
    string T;
    vector<int> dT;
    TlessDNA() : T("") {}
};

static TlessDNA MakeTless(const string& P)
{
    TlessDNA result;
    int dTaccum = 0;
    for(size_t i = 0; i < P.size(); i++)
    {
    	while(i < P.size() && P[i] == 'T') { dTaccum++; i++; }
    	result.dT.push_back(dTaccum);
	    dTaccum = 0;
    	if(i < P.size())
    	   result.T += P[i];
    }
    if(result.dT.size() == result.T.size())
	   result.dT.push_back(0);
    return result;
}

static string PrintTless(const TlessDNA& Tl)
{
    string R = ""; size_t i = 0;
    for(; i < Tl.T.size(); i++)
    {
        for(size_t t = 0; t < Tl.dT[i]; t++) R += "T";
        R += Tl.T[i];
    }
    if(i == Tl.T.size())
        for(size_t t = 0; t < Tl.dT[i]; t++) R += "T";

    return R;
}



#endif
