// by J@ :)
#ifndef OPTIONPARSER_H
#define OPTIONPARSER_H

#include "dna_util.h"
#include <iostream>

using namespace std;
using namespace libdna;


// defaults

const char* g_taligner_version = "v.3.0.1";
const string default_output_prefix = "TA3_output_";
const int default_num_threads = 10;

const int default_align_mm_max = 2;
const int default_align_seed_length = 8;
const int default_align_seed_step = 8;
const int default_align_max_indel_size = 16;
const double default_align_read_fraction_min = 0.75;
const int default_align_read_length_min = 14;
const int default_og_min_overlap = 10;
const int default_start_codon_offset = 0;  // left part
const int default_stop_codon_offset = 0; // right part
const double default_subsampling = 0; // off


struct TAlignerOptions
{
    vector<string> fastq_input_files;
    string output_prefix;
    bool print_svg;
    bool print_pdf;
    int num_threads;
    int job_size;

    int align_mm_max;
    int align_seed_length;
    int align_seed_step;
    int align_max_indel_size;

    double align_read_fraction_min;
    int align_read_length_min;

    int start_codon_offset;
    int stop_codon_offset;
    bool start_with_ATG;
    bool draw_orf_incomplete;
    int main_orf_id;
    int orf_length_min;

    double subsample_reads_fraction;
    int readgraph_min_overlap;
    vector<string> fasta_genes_files;

    bool dump_mapped_reads;

    bool draw_translatability_plot;

    vector<int> draw_these_orfs;
    vector<vector<int> > draw_orfs_color_rgb;

    bool btlc_draw;
    int btlc_pos;
    int btlc_k;
    int btlc_count;
};

void Parse_Arguments(int argc, char** argv, TAlignerOptions& options)
{
    int arg_num = 0;

    /*
        presetting defaults
    */

    options.print_svg = false;
    options.print_pdf = false;
    options.start_with_ATG = false;
    options.draw_orf_incomplete = true;
    options.main_orf_id = 0;
    options.orf_length_min = 0;
    options.num_threads = default_num_threads;
    options.output_prefix = default_output_prefix;

    options.align_mm_max = default_align_mm_max;
    options.align_seed_length = default_align_seed_length;
    options.align_seed_step = default_align_seed_step;
    options.align_max_indel_size = default_align_max_indel_size;
    options.align_read_fraction_min = default_align_read_fraction_min;
    options.align_read_length_min = default_align_read_length_min;
    options.start_codon_offset = default_start_codon_offset;
    options.stop_codon_offset = default_stop_codon_offset;
    options.subsample_reads_fraction = default_subsampling;
    options.readgraph_min_overlap = default_og_min_overlap;
    options.dump_mapped_reads = false;
    options.draw_translatability_plot = false;
    options.btlc_draw=false;
    options.btlc_pos=0;
    options.btlc_k=0;
    options.btlc_count=0;


    /*
        reading user's options
    */

    while(++arg_num < argc)
    {

        if(string(argv[arg_num]) == "--in_ref")
        {
			options.fasta_genes_files.push_back(string(argv[arg_num+1]));
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--in_lib")
        {
			options.fastq_input_files.push_back(string(argv[arg_num+1]));
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--out")
        {
			options.output_prefix = string(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--svg")
        {
			options.print_svg = true;
        }

        else if(string(argv[arg_num]) == "--pdf")
        {
			options.print_pdf = true;
        }

        else if(string(argv[arg_num]) == "--draw_transplot")
        {
            options.draw_translatability_plot = true;
        }

        else if(string(argv[arg_num]) == "--t" || string(argv[arg_num]) == "--num_threads")
        {
			options.num_threads = atoi(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--j" || string(argv[arg_num]) == "--job_size")
        {
			options.job_size = atoi(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--sl" || string(argv[arg_num]) == "--seed_length")
        {
			options.align_seed_length = atoi(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--ss" || string(argv[arg_num]) == "--seed_step")
        {
			options.align_seed_step = atoi(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--mm" || string(argv[arg_num]) == "--mismatch_max")
        {
			options.align_mm_max = atoi(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--mf" || string(argv[arg_num]) == "--min_fraction")
        {
			options.align_read_fraction_min = atof(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--mr" || string(argv[arg_num]) == "--min_read_size")
        {
			options.align_read_length_min = atoi(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--xi" || string(argv[arg_num]) == "--max_indel")
        {
			options.align_max_indel_size = atoi(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--start_codon_in")
        {
			options.start_codon_offset = atoi(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--stop_codon_in")
        {
			options.stop_codon_offset = atoi(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--start_with_atg")
        {
			options.start_with_ATG = true;
        }

        else if(string(argv[arg_num]) == "--subsample_reads")
        {
			options.subsample_reads_fraction = atof(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--rg_overlap_min")
        {
			options.readgraph_min_overlap = atoi(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--orf_length_min")
        {
            options.orf_length_min = atoi(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--report_orf_incomplete")
        {
			options.draw_orf_incomplete = false;
        }

        else if(string(argv[arg_num]) == "--set_main_orf")
        {
			options.main_orf_id = atoi(argv[arg_num+1]);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--dump_mapped_reads")
        {
			options.dump_mapped_reads = true;
        }

        else if(string(argv[arg_num]) == "--draw_orf")
        {
			options.draw_these_orfs.push_back(atoi(argv[arg_num+1]));
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--draw_orf_rgb")
        {
            auto params = Tokenize(string(argv[arg_num+1]),',');
			options.draw_these_orfs.push_back(atoi(params[0].c_str()));
            vector<int> color(3,0);
            color[0] = atoi(params[1].c_str());
            color[1] = atoi(params[2].c_str());
            color[2] = atoi(params[3].c_str());
            options.draw_orfs_color_rgb.push_back(color);
            arg_num++;
        }

        else if(string(argv[arg_num]) == "--draw_orfs_btlc")
        {
            auto params = Tokenize(string(argv[arg_num+1]), ',');

			options.btlc_pos = atoi(params[0].c_str());
            options.btlc_k = atoi(params[1].c_str());
            options.btlc_count = atoi(params[2].c_str());
            options.btlc_draw = true;
            arg_num++;
        }

    }

    if(options.draw_these_orfs.size() == 0)
        options.draw_these_orfs.push_back(0);
}


#endif
