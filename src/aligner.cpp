// by J@:)
#include "aligner.h"

#include <iostream>

vector<Subalignment> Do_Align_Read(ReadAlignmentTask T)
{
    vector<Subalignment> Alignments;

    // Critial modification!
    // was
    //for(size_t pos = 0; pos < T.read.T.size() - 2*T.seed; pos += T.seed_step)
    size_t SeedsPerRead = static_cast<size_t>(T.seed * (static_cast<int>(T.read.T.size()) / T.seed));
    for(size_t pos = 0; pos < SeedsPerRead - T.seed; pos += T.seed_step)
    {
        auto hits = Match_Index(T.read.T.substr(pos, T.seed), *(T.refIndex));

        for(auto& hit : hits)
        {
            bool appended = false;
            for(auto& align : Alignments)
            {
                if(hit.i == align.ref && hit.p == align.rd_rightmost + T.seed_step)
                {
                    appended = true;
                    align.rd_rightmost += T.seed_step;
                    break;
                }
            }

            if(!appended)
            {
                Subalignment align;
                align.rd_leftmost = pos;
                align.rd_rightmost = pos + T.seed_step;
                align.rf_start = hit.p;
                align.ref = hit.i;
                Alignments.push_back(align);
            }
        }
    }

    sort(Alignments.begin(), Alignments.end(), [&](Subalignment a, Subalignment b)
		{ return (a.rd_rightmost-a.rd_leftmost) > (b.rd_rightmost-b.rd_leftmost); });


    vector<Subalignment> Best;
    if(Alignments.size() > 0)
    {
        bool too_long_indel = false;
        for(size_t pos = 0; pos < T.read.dT.size(); pos++)
            if(T.read.dT[pos] > T.max_long_indel_size) {
                too_long_indel=true; break;
            }
        if(!too_long_indel)
            Best.push_back(Alignments[0]);
    }
    return Best;
}

vector<ReadAlignment> Filter_Invalid_Alignments(
    vector<Subalignment>& Alignments,
    ReadAlignmentTask& T)
{

    vector<ReadAlignment> ValidAlignments;

    for(size_t i = 0; i < Alignments.size(); i++)
    {
        ReadAlignment valid_alignment;
        valid_alignment.rf_start = Alignments[i].rf_start;
        valid_alignment.rd_from = Alignments[i].rd_leftmost;
        valid_alignment.rd_to = Alignments[i].rd_rightmost;
        valid_alignment.ref = Alignments[i].ref;

        const string& ref = (*T.Texts)[Alignments[i].ref].T;
        const string& read = T.read.T;

        // try to extend left
        int mm = 0;

        while(mm < T.mm &&
            valid_alignment.rf_start > 0 &&
            valid_alignment.rd_from > 0) {
            if(ref[--valid_alignment.rf_start] !=
            read[--valid_alignment.rd_from]) mm++;
        }

        int bonus_shift_left = 0;
        int bonus_shift_right = 0;

        while(mm < T.mm)
        {
            int pos_shift_left = 0;
            int pos_shift_right = 0;

            if(bonus_shift_left == 0)
            {
                while(valid_alignment.rf_start-pos_shift_left > 0 &&
                    valid_alignment.rd_from-pos_shift_left > 0 &&
                    ref[valid_alignment.rf_start-pos_shift_left] ==
                    read[valid_alignment.rd_from-pos_shift_left])
                    pos_shift_left++;
                bonus_shift_left = pos_shift_left;
            }

            if(bonus_shift_right == 0)
            {
                while(valid_alignment.rf_start+valid_alignment.rd_to+
                    pos_shift_right < ref.size() &&
                    valid_alignment.rd_to+pos_shift_right > 0 &&
                    ref[valid_alignment.rf_start+valid_alignment.rd_to+
                        pos_shift_right] ==
                    read[valid_alignment.rd_to+pos_shift_right])
                    pos_shift_right++;
                bonus_shift_right = pos_shift_right;
            }

            mm++;

            if(bonus_shift_left > bonus_shift_right)
            {
                valid_alignment.rf_start -= bonus_shift_left;
                valid_alignment.rd_from -= bonus_shift_left;
                bonus_shift_left = 0;
            }
            else
            {
                valid_alignment.rd_to += bonus_shift_right;
                bonus_shift_right = 0;
            }
        }

        if(T.align_rd_min*read.size() < valid_alignment.rd_to-valid_alignment.rd_from)
        {

            valid_alignment.read = T.read;
            ValidAlignments.push_back(valid_alignment);
        }

    }
    return ValidAlignments;
}


bool Read_Alignment_Matches_Reference(
    ReadAlignment& Alignment,
    vector<TlessDNA>& Texts)
{
        int unedited_positions = 0;
        for(size_t i = 1; i < Alignment.rd_to - Alignment.rd_from; i++)
        {
            if(Alignment.read.dT[Alignment.rd_from+i] !=
                Texts[Alignment.ref].dT[Alignment.rf_start+i])
                {
                    return false;
                }
        }
        Alignment.as_ref = true;
        return true;

}


vector<AlignedRow> Calculate_Alignment_Matrix(
    vector<ReadAlignment>& Alignments,
    const TlessDNA& reference)
{
    sort(Alignments.begin(), Alignments.end(), [&](ReadAlignment a, ReadAlignment b)
		{ return a.rf_start < b.rf_start; });

    vector<AlignedRow> matrix;
    for(auto& alignment : Alignments)
    {
        AlignedRow row;
        row.ref = alignment.ref;
        row.rf_start = alignment.rf_start;
        row.read = &alignment.read;
        row.rd_from = alignment.rd_from;
        for(size_t i = 0; i < alignment.rd_to-alignment.rd_from; i++)
        {
            row.r.push_back(alignment.read.dT[alignment.rd_from+i+1] -
                reference.dT[alignment.rf_start+i+1]);
        }
        matrix.push_back(row);
    }
    return matrix;
}

vector<MappedPart> Truncate_Unmapped_Parts(vector<ReadAlignment>& Alignments)
{
    vector<MappedPart> mapped_parts;
    for(auto& align : Alignments)
    {
        MappedPart MP;
        MP.as_ref = align.as_ref;
        MP.rf_start = align.rf_start;
        MP.ref = align.ref;
        MP.supp = 0;
        MP.mp.T = align.read.T.substr(align.rd_from, align.rd_to-align.rd_from);
        for(int i = 0; i <= align.rd_to-align.rd_from; i++)
        {
            MP.mp.dT.push_back(align.read.dT[align.rd_from+i]);
        }
        mapped_parts.push_back(MP);
    }
    return mapped_parts;
}
