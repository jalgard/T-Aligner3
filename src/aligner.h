// by J@:)

#ifndef ALIGNER_H
#define ALIGNER_H

#include "indexref.h"
#include "tless.h"
#include "dna_util.h"

using namespace libdna;
using namespace std;


struct ReadAlignment
{
    size_t rf_start;
    size_t ref;
    size_t rd_from;
    size_t rd_to;

    TlessDNA read;
    bool as_ref;
};

struct MappedPart
{
    TlessDNA mp;
    size_t rf_start;
    size_t ref;
    int supp;
    bool as_ref;
};

struct Subalignment
{
    size_t rd_leftmost;
    size_t rd_rightmost;
    size_t rf_start;
    size_t ref;
};

struct ReadAlignmentTask
{
    unordered_map<string, vector<Rindex> >* refIndex;
    vector<TlessDNA>* Texts;
    TlessDNA read;
    int seed; // seed length
    int seed_step; // seed step
    int mm; // total possible number of mismatches
    int max_long_indel_size; // maximal number of T's in read
    double align_rd_min; // minimal fraction of read to be aligned
    int min_rd_length;
};

vector<Subalignment> Do_Align_Read(ReadAlignmentTask T);

vector<ReadAlignment> Filter_Invalid_Alignments(
    vector<Subalignment>& Alignments,
    ReadAlignmentTask& T);

bool Read_Alignment_Matches_Reference(
    ReadAlignment& Alignment,
    vector<TlessDNA>& Texts);


struct AlignedRow
{
    vector<int> r;
    size_t ref;
    size_t rf_start;
    TlessDNA* read;
    size_t rd_from;
};

vector<AlignedRow> Calculate_Alignment_Matrix(
    vector<ReadAlignment>& Alignments,
    const TlessDNA& reference);

vector<MappedPart> Truncate_Unmapped_Parts(vector<ReadAlignment>& Alignments);

#endif
