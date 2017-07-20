// by J@ :)


#ifndef ORFFINDER_H
#define ORFFINDER_H

#include "dna_util.h"
#include "option_parser.h"
#include "aligner.h"

#include <set>

using namespace std;
using namespace libdna;

/*
    vector<ReadAlignment> Alignments
    vector<TlessDNA> Refbase

    struct ReadAlignment:
    (size_t) rf_start, ref, rd_from, rd_to;
    (TlessDNA) read;
    (bool) as_ref;
*/


void Assmeble_Unedited_Alignments(
    vector<MappedPart>& Alignments,
    const TlessDNA& Reference,
    vector<MappedPart>& EditedAlignments,
    vector<int>& Reference_coverage
)
{
    Reference_coverage.resize(Reference.T.size(),0);

    for(auto& alignment : Alignments)
    {
        if(alignment.as_ref)
        {
            for(int p = 0; p < alignment.mp.T.size(); p++)
            {
                Reference_coverage[alignment.rf_start+p]++;
            }
        }
        else
        {
            EditedAlignments.push_back(alignment);
        }
    }
}


vector<MappedPart> Join_Exact_Mappings(vector<MappedPart>& mapped)
{
    map<size_t, int> skip;
    vector<MappedPart> picked;

    for(size_t i = 0; i < mapped.size(); i++)
    {
        if(skip.count(i)) continue;
        int supp_counter = 0;

        for(size_t j = i+1; j < mapped.size(); j++)
        {
            if(mapped[i].rf_start == mapped[j].rf_start &&
                mapped[i].mp.T.size() == mapped[j].mp.T.size() &&
                mapped[i].mp.dT.size() == mapped[j].mp.dT.size())
            {
                bool equal = true;
                for(int k = 0; k < mapped[i].mp.dT.size(); k++)
                {
                    if(mapped[i].mp.dT[k] != mapped[j].mp.dT[k])
                    {
                        equal = false;
                        break;
                    }
                }
                if(equal)
                {
                    skip[j] = 1;
                    supp_counter++;
                }
            }
        }

        mapped[i].supp = supp_counter+1;
        picked.push_back(mapped[i]);
    }
    return picked;
}


vector<MappedPart> Join_Substrings(vector<MappedPart>& mapped)
{
    map<size_t, int> skip;
    vector<MappedPart> picked;

    // sort by refstart and then by size (rev order)

    stable_sort(mapped.begin(), mapped.end(), [](MappedPart a, MappedPart b)
        { return a.rf_start < b.rf_start; });
    stable_sort(mapped.begin(), mapped.end(), [](MappedPart a, MappedPart b)
    	{ return a.mp.T.size() > b.mp.T.size(); });
    for(size_t i = 0; i < mapped.size(); i++)
    {
        if(skip.count(i)) continue;
        int supp_counter = 0;
        for(size_t j = i+1; j < mapped.size(); j++)
        {
            if(mapped[j].rf_start >= mapped[i].rf_start &&
            mapped[j].rf_start+mapped[j].mp.T.size() <=
            mapped[i].rf_start+mapped[i].mp.T.size())
            {
                bool is_substring = true;
                int j_i_shift = mapped[j].rf_start-mapped[i].rf_start;
                for(int k = 0; k < mapped[j].mp.dT.size(); k++)
                {
                    if(mapped[i].mp.dT[j_i_shift+k] != mapped[j].mp.dT[k])
                    {
                        is_substring = false;
                        break;
                    }
                }
                if(is_substring)
                {
                    skip[j] = 1;
                    supp_counter++;
                }
            }
        }

        mapped[i].supp = supp_counter+1;
        picked.push_back(mapped[i]);
    }

    return picked;
}


static MappedPart Join_Mapped_Parts(MappedPart L, MappedPart R)
{
    MappedPart M; M.supp = -1;
    if(R.rf_start >= L.rf_start+L.mp.T.size()) return M;
    if(R.rf_start < L.rf_start) return M;
    int l_r_shift = R.rf_start - L.rf_start;
    int overlap_length = L.mp.T.size()-l_r_shift;
    bool join = true;
    for(int k = 0; k < overlap_length; k++)
    {
        if(L.mp.dT[l_r_shift+k] != R.mp.dT[k])
        {
            join = false;
            break;
        }
    }
    if(join)
    {
        M.mp.T = L.mp.T.substr(0, l_r_shift) + R.mp.T;
        for(int i = 0; i < l_r_shift; i++)
            M.mp.dT.push_back(L.mp.dT[i]);
        for(int i = 0; i < R.mp.dT.size(); i ++)
            M.mp.dT.push_back(R.mp.dT[i]);
        M.rf_start = L.rf_start;
        M.ref = L.ref;
        M.supp = L.supp+R.supp;
    }

    return M;
}

struct OverlapEdge
{
    int ddist;
    int8_t tdist;
};


// OG (overlap graph) is valid ONLY for returned MPN (mapped part nodes) order
void Build_Overlap_Graph(
    vector<MappedPart>& MPN,
    map<int, map<int, OverlapEdge> >& OG,
    int min_overlap = 10,
    int min_dadd = 10)
{
    cout << "Building overlap graph...\n";
    OG.clear();
    sort(MPN.begin(), MPN.end(), [](MappedPart a, MappedPart b)
        { return a.mp.T.size() > b.mp.T.size(); });
    stable_sort(MPN.begin(), MPN.end(), [](MappedPart a, MappedPart b)
		{ return a.rf_start < b.rf_start; });


    int edges = 0;
    for(int i = 0; i < MPN.size(); i++)
    {
        for(int j = 0; j < MPN.size(); j++)
        {
            if(i==j) continue;
            if(MPN[j].rf_start + min_overlap <
                MPN[i].rf_start + MPN[i].mp.T.size())
            {
                int l_r_shift = MPN[j].rf_start - MPN[i].rf_start;
                int overlap_length = MPN[i].mp.T.size()-l_r_shift;
                bool join = true;
                for(int k = 0; k < overlap_length; k++)
                {
                    if(l_r_shift + k < MPN[i].mp.dT.size() &&
                        k < MPN[j].mp.dT.size() &&
                        MPN[i].mp.dT[l_r_shift+k] != MPN[j].mp.dT[k])
                    {
                        join = false;
                        break;
                    }
                }
                if(join)
                {
                    OverlapEdge E; E.tdist = 0; E.ddist = 0;
                    OG[i][j] = E;
                    edges++;
                }
            }
        }
    }
    cout << "Overlap graph done with " << edges << " edges\n";

}

vector<int> Select_Reads_With_Start_Codon(
    vector<MappedPart>& MPN,
    const TlessDNA& Reference,
    int offset=0,
    bool start_w_atg = false)
{
    set<int> with_start_nodes;
    vector<int> result;
    for(int n = 0; n < MPN.size(); n++)
    {
        if(offset > 0 && MPN[n].rf_start > offset) continue;
        else if(MPN[n].rf_start > Reference.T.size()/2) continue;
        string read = PrintTless(MPN[n].mp);
        for(int f = 0; f < 3; f++)
        {
            bool has_start = false;
            for(int i = 0; i < read.size()-6; i += 3)
            {
                string codon = read.substr(f+i, 3);
                if(codon == "ATG" || (codon == "ATA" && !start_w_atg)) has_start = true;
                if(codon == "TAG" || codon == "TAA") has_start = false;
            }
            if(has_start) with_start_nodes.insert(n);
        }
    }
    for(auto& x : with_start_nodes) result.push_back(x);
    sort(result.begin(), result.end(), [&](int a, int b)
        {return MPN[a].rf_start < MPN[b].rf_start; });
    return result;
}

vector<int> Select_Reads_With_Stop_Codon(
    vector<MappedPart>& MPN,
    const TlessDNA& Reference,
    int offset = 0)
{
    set<int> with_stop_nodes;
    vector<int> result;
    for(int n = 0; n < MPN.size(); n++)
    {
        if(offset > 0 && MPN[n].rf_start < Reference.T.size()-offset) continue;
        else if(MPN[n].rf_start < Reference.T.size()/2) continue;
        string read = PrintTless(MPN[n].mp);
        for(int f = 0; f < 3; f++)
        {
            bool has_stop = false;
            for(int i = 0; i < read.size()-6; i += 3)
            {
                string codon = read.substr(f+i, 3);
                if(codon == "TAG" || codon == "TAA") has_stop = true;
            }
            if(has_stop) with_stop_nodes.insert(n);
        }
    }
    for(auto& x : with_stop_nodes) result.push_back(x);
    sort(result.begin(), result.end(), [&](int a, int b)
        {return MPN[a].rf_start < MPN[b].rf_start; });
    return result;
}

static string Find_Longest_ORF(MappedPart& M, bool start_w_atg=false)
{
    string dnaseq = PrintTless(M.mp);

    int lstart = 0; int lend = 0;
    for(int f = 0; f < 3; f++)
    {
        int cstart = 0; int cend = 0;
        for(int i = 0; i < dnaseq.size() - 6; i += 3)
        {
            string codon = dnaseq.substr(f+i, 3);
            if(codon == "ATG" || (codon == "ATA" && !start_w_atg) ) {
                if(cstart == 0)
                    cstart = f+i;
            }
            if(codon == "TAG" || codon == "TAA") {
                cend = f+i;
                if(cstart > 0 && (cend - cstart) > (lend - lstart)) {
                    lend = cend;
                    lstart = cstart;
                }
                cstart = 0;
                cend = 0;
            }
        }
        if(cstart > 0 && cend == 0 && (dnaseq.size() - cstart) > (lend - lstart) )
        {
            lend = dnaseq.size(); lstart = cstart;
        }
    }

    return dnaseq.substr(lstart, lend-lstart);

}

static string Find_Longest_ORF_Complete(MappedPart& M, bool start_w_atg=false)
{
    string dnaseq = PrintTless(M.mp);

    int lstart = 0; int lend = 0;
    for(int f = 0; f < 3; f++)
    {
        int cstart = 0; int cend = 0;
        for(int i = 0; i < dnaseq.size() - 6; i += 3)
        {
            string codon = dnaseq.substr(f+i, 3);
            if(codon == "ATG" || (codon == "ATA" && !start_w_atg)) {
                if(cstart == 0)
                    cstart = f+i;
            }
            if(codon == "TAG" || codon == "TAA") {
                cend = f+i;
                if(cstart > 0 && (cend - cstart) > (lend - lstart)) {
                    lend = cend;
                    lstart = cstart;
                }
                cstart = 0;
                cend = 0;
            }
        }
    }

    return dnaseq.substr(lstart, lend-lstart);

}

struct taORF
{
    string orf;
    int tlstart;
    int tlend;
};

struct floParam
{
    MappedPart M;
    bool start_w_atg;
    bool is_complete;
};

static taORF Find_Longest_ORF(floParam flp)
{
    string dnaseq = PrintTless(flp.M.mp);

    int lstart = 0; int lend = 0;
    for(int f = 0; f < 3; f++)
    {
        int cstart = 0; int cend = 0;
        for(int i = 0; i < dnaseq.size() - 6; i += 3)
        {
            string codon = dnaseq.substr(f+i, 3);
            if(codon == "ATG" || (codon == "ATA" && flp.start_w_atg==false)) {
                if(cstart == 0)
                    cstart = f+i;
            }
            if(codon == "TAG" || codon == "TAA") {
                cend = f+i;
                if(cstart > 0 && (cend - cstart) > (lend - lstart)) {
                    lend = cend;
                    lstart = cstart;
                }
                cstart = 0;
                cend = 0;
            }
        }
        if(flp.is_complete==false && cstart > 0 && cend == 0 &&
            (dnaseq.size() - cstart) > (lend - lstart))
        {
            lend = dnaseq.size(); lstart = cstart; // mode may2017: "lend = dnaseq.size()-1" influenced ORF order
        }
    }

    taORF tad;
    tad.orf = dnaseq.substr(lstart, lend-lstart);

    int t_ = 0;
    for(int i = 0; i <= lstart; i++)
        if(dnaseq[i] != 'T') t_++;
    tad.tlstart = t_;
    t_=0;
    for(int i = 0; i <= lend; i++)
        if(dnaseq[i] != 'T') t_++;
    tad.tlend = t_;

    return tad;
}


MappedPart BFS_For_LongestORF(
    int source,
    vector<MappedPart>& MPN,
    map<int, map<int, OverlapEdge> >& OG,
    vector<int> Starts,
    vector<int> Stops,
    int min_overlap = 10,
    bool start_with_ATG = false
)
{
    MappedPart M; M.supp = -1;
    vector<int> dist(MPN.size(), 0);
    vector<MappedPart> prev(MPN.size(), M);
    vector<int> Q; Q.push_back(source);
    vector<int> seen(MPN.size(),0);

    prev[source] = MPN[source];

    while(Q.size() > 0)
    {
        int node = Q.back();
        Q.pop_back();

        MappedPart MaxPossibleExt;
        int MaxExtDist = 0;
        int MaxExtNode = 0;

        for(auto& nb : OG[node])
        {
            //if(seen[nb.first] > 0) continue;
            int CurrentExtDist = 0; MappedPart CurrentPossibleExt;
            if(prev[node].rf_start +  min_overlap < MPN[nb.first].rf_start &&
                MPN[nb.first].rf_start + min_overlap < prev[node].rf_start +
                prev[node].mp.T.size())
            {
                CurrentPossibleExt = Join_Mapped_Parts(prev[node], MPN[nb.first]);
                if(CurrentPossibleExt.mp.T.size() > min_overlap)
                {
                    floParam flp;
                    flp.M = CurrentPossibleExt;
                    flp.start_w_atg = start_with_ATG;
                    flp.is_complete = false;
                    //if(OG.count(nb.first)==0) flp.is_complete = true;
                    CurrentExtDist = Find_Longest_ORF(flp).orf.size();
                }
            }
            if(CurrentExtDist > MaxExtDist ||
            (CurrentExtDist == MaxExtDist &&
            find(Stops.begin(), Stops.end(), nb.first) != Stops.end()))
            {
                MaxExtDist = CurrentExtDist;
                MaxExtNode = nb.first;
                MaxPossibleExt = CurrentPossibleExt;
            }
        }

        if(dist[MaxExtNode] < MaxExtDist)
        {
            dist[MaxExtNode] = MaxExtDist;
            //seen[MaxExtNode] = 1;
            prev[MaxExtNode] = MaxPossibleExt;
            Q.insert(Q.begin(), MaxExtNode);
        }
    }

    int max_orf = 0; int max_orf_id = 0;
    for(auto& s : Stops)
    {
        if(max_orf < dist[s]) {
            max_orf = dist[s];
            max_orf_id = s;
        }
    }

    return prev[max_orf_id];

}

static string TranslateORF(const std::string& sequence)
{
    string translation = "";
    map<string, string> code1;
	code1["TTT"]="F";	code1["TCT"]="S";	code1["TAT"]="Y";	code1["TGT"]="C";
	code1["TTC"]="F";	code1["TCC"]="S";	code1["TAC"]="Y";	code1["TGC"]="C";
	code1["TTA"]="L";	code1["TCA"]="S";	code1["TAA"]="*";	code1["TGA"]="W";
	code1["TTG"]="L";	code1["TCG"]="S";	code1["TAG"]="*";	code1["TGG"]="W";
	code1["CTT"]="L";	code1["CCT"]="P";	code1["CAT"]="H";	code1["CGT"]="R";
	code1["CTC"]="L";	code1["CCC"]="P";	code1["CAC"]="H";	code1["CGC"]="R";
	code1["CTA"]="L";	code1["CCA"]="P";	code1["CAA"]="Q";	code1["CGA"]="R";
	code1["CTG"]="L";	code1["CCG"]="P";	code1["CAG"]="Q";	code1["CGG"]="R";
	code1["ATT"]="I";	code1["ACT"]="T";	code1["AAT"]="N";	code1["AGT"]="S";
	code1["ATC"]="I";	code1["ACC"]="T";	code1["AAC"]="N";	code1["AGC"]="S";
	code1["ATA"]="I";	code1["ACA"]="T";	code1["AAA"]="K";	code1["AGA"]="R";
	code1["ATG"]="M";	code1["ACG"]="T";	code1["AAG"]="K";	code1["AGG"]="R";
	code1["GTT"]="V";	code1["GCT"]="A";	code1["GAT"]="D";	code1["GGT"]="G";
	code1["GTC"]="V";	code1["GCC"]="A";	code1["GAC"]="D";	code1["GGC"]="G";
	code1["GTA"]="V";	code1["GCA"]="A";	code1["GAA"]="E";	code1["GGA"]="G";
	code1["GTG"]="V";	code1["GCG"]="A";	code1["GAG"]="E";	code1["GGG"]="G";
    for(int i = 0; i < sequence.size(); i+=3)
        translation += code1[sequence.substr(i,3)];

    return translation;
}






#endif
