// by J@:)
#ifndef INDEXREF_H
#define INDEXREF_H

#include "tless.h"
#include <unordered_map>


using namespace std;

struct Rindex
{
    size_t p;
    size_t i;
};

// Build hash from collection of Texts
static void Build_Fasta_Index(vector<TlessDNA>& Texts, unordered_map<string, vector<Rindex> >& Index, int seed = 8)
{
    for(size_t T = 0; T < Texts.size(); T++)
    {
        for(size_t p = 0; p < Texts[T].T.size() - seed; p++) // extract seed with minimal possible step 1
        {
            Rindex ri; ri.p = p; ri.i = T;
            Index[Texts[T].T.substr(p, seed)].push_back(ri);
        }
    }
}

// Look for a Seed match in Index
// S should be size of seed (8)
static vector<Rindex> Match_Index(const string& S, unordered_map<string, vector<Rindex> >& Index)
{
    auto hit = Index.find(S);
    if(hit != Index.end())
    {
        return hit->second;
    }
    vector<Rindex> M;
    return M;
}

#endif
