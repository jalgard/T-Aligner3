// by J@, 2016

#ifndef CALCEDITING_H
#define CALCEDITING_H

#include <future>
#include "aligner.h"
#include "option_parser.h"

void Calculate_Editing_Pyramid(vector<AlignedRow>& alignedReads,
  map<int, int>& editingSitePyramid,
  map<int, int>& editingDistPyramid,
  map<string, int>& statsOut, int mappedSegment = 100)
{

    int skipped_too_short = 0;
    for(auto& aRead: alignedReads)
    {
        if(aRead.r.size() < mappedSegment)
        {
            skipped_too_short++; continue;
        }

        int edited_pos_count = 0;
        int edited_dist = 0;
        for(int i = 0; i < mappedSegment; i++)
        {
            if(aRead.r[i] != 0)
            {
                edited_pos_count++;
                edited_dist += abs(aRead.r[i]);
            }

        }
        editingSitePyramid[edited_pos_count]++;
        editingDistPyramid[edited_dist]++;

    }
    statsOut["SkippedTooShort"] = skipped_too_short;
}


void Find_New_Cryptogenes(vector<AlignedRow>& alignedReads,
    map<int, int>& editedSegments, int binSize=200)
{
    for(auto& aRead: alignedReads)
    {
        int epc = 0;
        for(int i = 0; i < aRead.r.size(); i++)
        {
            if(aRead.r[i] != 0) epc++;
        }
        if(epc > 10)
        {
            editedSegments[aRead.rf_start / 200]++;
        }
    }
}

#endif
