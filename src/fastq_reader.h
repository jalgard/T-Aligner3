// by J@, 2016

#ifndef FASTQREADER_H
#define FASTQREADER_H

#include<future>
#include "aligner.h"
#include "option_parser.h"

inline void Parse_One_Alignment(
    vector<ReadAlignment>& Alignments,
    ReadAlignmentTask Task)
{
    if(Task.read.T.size() < Task.min_rd_length) return;
    vector<Subalignment> seedAlignments = Do_Align_Read(Task);

    vector<ReadAlignment> finalAlignments =
        Filter_Invalid_Alignments(seedAlignments, Task);
    if(finalAlignments.size() > 0)
        Alignments.push_back(finalAlignments[0]);
}

inline vector<ReadAlignment> Do_Alignments_Chain(vector<ReadAlignmentTask> JobChain)
{
    vector<ReadAlignment> Alignments;
    for(size_t i = 0; i < JobChain.size(); i++)
    {
        Parse_One_Alignment(Alignments, JobChain[i]);
    }
    return Alignments;
}


//
vector<ReadAlignment> Align_Reads_From_Fastq(
    TAlignerOptions Opts, int lib_number,
    ReadAlignmentTask& AlignmentOptions)
{

    string l1,l2,l3,l4;
    ifstream ifile(Opts.fastq_input_files[lib_number].c_str());
    vector<ReadAlignment> Alignments;

    while(1)
    {
        vector<string> reads_read;
        while(
            getline(ifile, l1, '\n') &&
            getline(ifile, l2, '\n') &&
            getline(ifile, l3, '\n') &&
            getline(ifile, l4, '\n'))
        {
            reads_read.push_back(l2);
            if(reads_read.size()==Opts.num_threads*Opts.job_size) break;
        }
        if(reads_read.size() < Opts.num_threads) break;

        // at least num_threads reads we have - do fow/rev alignment
        vector<future<vector<ReadAlignment> > > align_results;
        for(int rd = 0; rd < reads_read.size(); rd += Opts.job_size)
        {
            vector<ReadAlignmentTask> CurrentBatchJob;
            for(int j = 0; j < Opts.job_size; j++)
            {
                if(rd + j >= reads_read.size()) break;
                ReadAlignmentTask Task_fow = AlignmentOptions;
                Task_fow.read = MakeTless(reads_read[rd+j]);
                if(Task_fow.read.T.size() > Opts.align_seed_length * 3)
                    CurrentBatchJob.push_back(Task_fow);
                ReadAlignmentTask Task_rev = Task_fow;
                Task_rev.read = MakeTless(rcDNA(reads_read[rd+j]));
                if(Task_rev.read.T.size() > Opts.align_seed_length * 3)
                    CurrentBatchJob.push_back(Task_rev);

            }
            align_results.push_back(async(launch::async,
                Do_Alignments_Chain, CurrentBatchJob));
        }

        for(auto &thread_align : align_results)
        {
        	auto result = thread_align.get();
        	Alignments.insert(Alignments.end(), result.begin(), result.end());
        }

    }
    return Alignments;
}

#endif
