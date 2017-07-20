// by J@:)

#ifndef ALIGNWRITER_H
#define ALIGNWRITER_H

#include "aligner.h"

#include <fstream>

void Do_Save_All_Reads_Alignment_FASTA(const char* file,
	const TlessDNA& reference,
	vector<ReadAlignment>& good_reads)
{

	vector<int> max_t(reference.T.size()+1,0);

	// scan ref for max dT
	for(int i=0; i < reference.T.size(); i++)
	{
		if(reference.dT[i+1] > max_t[i+1]) max_t[i+1] = reference.dT[i+1];
	}
	// scan reads for max dT //
	for(auto x : good_reads)
	{
        for(int i = 0; i < x.rd_to-x.rd_from; i++)
        {
			if(x.read.dT[x.rd_from+i] > max_t[x.rf_start+i])
                max_t[x.rf_start+i] = x.read.dT[x.rd_from+i];
		}
	}

	// write the alignment //
	ofstream ofile(file);

	// write ref //
    ofile << ">Reference\n";
	for(int i=0; i < reference.T.size(); i++)
	{
		ofile << reference.T[i];
		for(int j=0; j < reference.dT[i+1]; j++)
			ofile << "T";
		for(int j=reference.dT[i+1]; j < max_t[i+1]; j++)
			ofile << "-";
	}
	ofile << "\n";
	// write reads //
    stable_sort(good_reads.begin(), good_reads.end(), [&](ReadAlignment a, ReadAlignment b)
		{ return a.rf_start > b.rf_start; });

    int rd_counter = 0;
	for(auto& x : good_reads)
	{
		ofile << ">Read_" << (++rd_counter) << "_at_" << x.rf_start << "\n";

        for(int i=0; i < x.rf_start; i++)
		{
			ofile << "-";
			for(int j=0; j < max_t[i+1]; j++) ofile << "-";
		}

		for(int i=0; i < x.rd_to-x.rd_from-1 && i < x.read.dT.size()-x.rd_from-1; i++)
		{
			ofile << x.read.T[x.rd_from+i];
			for(int j=0; j < x.read.dT[x.rd_from+i+1]; j++)
				ofile << "T";
			for(int j=x.read.dT[x.rd_from+i+1]; j < max_t[x.rf_start+i+1]; j++)
				ofile << "-";
		}
		ofile << "\n";
	}

}
static void DumpTless(const TlessDNA& Tl, const char* to_file)
{
    ofstream ofile(to_file);

	int gpos = 0;
	for(int i = 0; i < Tl.dT.size(); i++)
	{
		int pos = i+1;
		string Ts = ""; for(int j = 0; j < Tl.dT[i]; j++) Ts += "T";
		char tll = Tl.T[i];
		gpos += 1; gpos += Ts.size();
		ofile << pos << "\t" << Ts << "\t" << tll << "\t" << gpos << "\n";
	}

	ofile.close();
}

#endif
