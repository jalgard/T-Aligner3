// by J@:)
#include "dna_util.h"
#include "tless.h"
#include "indexref.h"
#include "fastq_reader.h"
#include "option_parser.h"
#include "align_writer.h"
#include "dotmatrix_drawer.h"
#include "frames_drawer.h"
#include "orf_finder.h"

#include <QtGui>
//#include <QSvgGenerator>


using namespace std;
using namespace libdna;

int main(int argc, char** argv)
{
    cout << "Starting T-Aligner " << g_taligner_version << "\n";

    TAlignerOptions Program_Options;
    Parse_Arguments(argc, argv, Program_Options);

    QApplication app(argc, argv);
	QWidget* drawerWidget = new QWidget(0);
    QGraphicsView* graphView = new QGraphicsView(drawerWidget);
    QGraphicsScene* Scene = new QGraphicsScene();
    graphView->setScene(Scene);

    QGraphicsScene* MappedReadsScene = new QGraphicsScene();


    // Reference loading
    vector<vector<string> > Reference_data;
    vector<TlessDNA > Reference_Tless;
    for(auto& ref_input : Program_Options.fasta_genes_files)
        read_fasta(ref_input.c_str(), Reference_data);

    for(auto& x : Reference_data)
    {
        auto tl = MakeTless(x[1]);
        Reference_Tless.push_back(tl);
    }

    unordered_map<string, vector<Rindex> > Index;
    Build_Fasta_Index(Reference_Tless, Index, Program_Options.align_seed_length);

    ReadAlignmentTask AlignmentOptions;
    AlignmentOptions.refIndex = &Index;
    AlignmentOptions.seed = Program_Options.align_seed_length;
    AlignmentOptions.seed_step = Program_Options.align_seed_step;
    AlignmentOptions.mm = Program_Options.align_mm_max;
    AlignmentOptions.max_long_indel_size = Program_Options.align_max_indel_size;
    AlignmentOptions.align_rd_min = Program_Options.align_read_fraction_min;
    AlignmentOptions.min_rd_length = Program_Options.align_read_length_min;
    AlignmentOptions.Texts = &Reference_Tless;

    vector<ReadAlignment> Alll =
        Align_Reads_From_Fastq(Program_Options, 0, AlignmentOptions);


    // Align second fastq library
    vector<ReadAlignment> ALib2;
    if(Program_Options.mode_lib_compare)
    {
        ALib2 = Align_Reads_From_Fastq(Program_Options, 1, AlignmentOptions);
    }
    //cout << "Alignment done!\n\n";
    //Do_Save_All_Reads_Alignment_FASTA((Program_Options.output_prefix+"_Alignment.fasta").c_str(),
    //	Reference_Tless[0], Alll);

    if(Program_Options.dump_mapped_reads == true)
    {
        ofstream mapped_part_fastq((Program_Options.output_prefix+"_mapped.fastq").c_str());
        int mapped_read_counter = 0;
        for(int i = 0; i < Alll.size(); i++)
        {
            string nucs = PrintTless(Alll[i].read);
            mapped_part_fastq << "@rd" << ++mapped_read_counter << "\n" <<
                nucs << "\n+\n" << string(nucs.size(),'I') << "\n";
        }
    }

    int mmmcount = 0;

    vector<AlignedRow> matrix = Calculate_Alignment_Matrix(Alll,
        Reference_Tless[0]);

    // alignment matrix for lib 2
    vector<AlignedRow> matrix_lib2;
    if(Program_Options.mode_lib_compare)
    {
        matrix_lib2 = Calculate_Alignment_Matrix(ALib2, Reference_Tless[0]);
    }

    int rd = 0;
    for(auto& x : Alll)
    {
        if(Read_Alignment_Matches_Reference(
           x, Reference_Tless))
           {
               mmmcount++;
           }
           rd++;
    }


    cout << "Match ref\t" << mmmcount << endl;
    cout << "Aligned\t" << matrix.size() << endl;

    return 0;
    // EARLY EXIT

    QPointF upperleft_point(0,100);

    SceneFramesParam SFP;
    SFP.scene = Scene;
    SFP.UpperCorner = upperleft_point;
    SFP.n_frames = 5;
    SFP.x_axis_step = 7.0;
    SFP.frame_heights.push_back(200);
    SFP.frame_heights.push_back(200);
    SFP.frame_heights.push_back(200);
    SFP.frame_heights.push_back(13*(Program_Options.draw_these_orfs.size()-1));
    SFP.frame_heights.push_back(200);


    /*
    A, editing states
    B, coverage, editing rate
    C, translatable editing states, ORFs
    D, cloud coverage, ORFs
    */
    SFP.frame_names.push_back("A, editing states");
    SFP.frame_names.push_back("B, coverage,\n% edited reads");
    SFP.frame_names.push_back("C, translatable\nediting states, ORFs");
    SFP.frame_names.push_back("");
    SFP.frame_names.push_back("D, cloud coverage,\nORFs");


    if(Program_Options.btlc_draw || Program_Options.mode_lib_compare) {
        SFP.n_frames += 1;
        SFP.frame_heights.push_back(200);
        SFP.frame_names.push_back("E, selected ORFs");

    }

    //SFP.frame_heights.push_back(100);
    //SFP.frame_names.push_back("Edited reads\n mapped on ORFs");

    SFP.Reference = Reference_Tless[0];
    SFP.mapping_stats = "Reads mapped " + i2s(matrix.size()) + "\n";
    SFP.heading_label = "TAligner v.3.0";
    SFP.pixmap = QPixmap (":/logo.png");

    AbstractDrawOptions ADO_dots, ADO_cloud, ADO_orf, ADO_orf_tp, ADO_orf_cloud, ADO_cover;
    ADO_dots.Scene = Scene; ADO_dots.UpperCorner = SFP.GetFrameCenter(1);
    ADO_cover.Scene = Scene; ADO_cover.UpperCorner = SFP.GetFrameCenter(2);
    ADO_cloud.Scene = Scene; ADO_cloud.UpperCorner = SFP.GetFrameCenter(3) + QPointF(0, 40.0);
    //ADO_cloud.mColor = QColor(0,0,0,25);


    cout << "Total alignments gone to trimming step: " << Alll.size() << "\n";
    auto All = Truncate_Unmapped_Parts(Alll);
    Alll.clear();


    // Filters round
    cout << "Before exact filter " << All.size() << "\n";
    auto A1 = Join_Exact_Mappings(All); //All.clear();
    cout << "After exact filter " << A1.size() << "\n";
    auto A2 = Join_Substrings(A1); A1.clear();
    cout << "After substring filter " << A2.size() << "\n";

    vector<MappedPart> MPN = A2; A2.clear();



    if(Program_Options.subsample_reads_fraction > 0)
    {
        std::random_shuffle ( MPN.begin(), MPN.end() );
        MPN.resize(Program_Options.subsample_reads_fraction*MPN.size());
    }


    cout << "Reads selected for ORF assembly: " << MPN.size() << "\n";
    map<int, map<int, OverlapEdge> > OG;

    Build_Overlap_Graph(MPN, OG, Program_Options.readgraph_min_overlap,
        Program_Options.start_with_ATG);

    auto starts = Select_Reads_With_Start_Codon(MPN, Reference_Tless[0],
        Program_Options.start_codon_offset, Program_Options.start_with_ATG);
    cout << "Starts: " << starts.size() << "\n";

    auto stops = Select_Reads_With_Stop_Codon(MPN, Reference_Tless[0],
        Program_Options.stop_codon_offset);
    cout << "Stops: " << stops.size() << "\n";

    cout << "Done!\nStarting assembly...\n";

    vector<MappedPart> DFS_results;
    for(auto& s : starts)
        DFS_results.push_back(BFS_For_LongestORF(
            s, MPN, OG, starts, stops, 10, Program_Options.start_with_ATG));


    struct Orfhit
    {
        string aa;
        string nuc;
        int start, end;
        int supp;
        size_t id;
        int edit_dist_2_main;
        int total_edit_dist_2_main;
        int es_supp;
        double es_ratio;
        int mrna_edit_dist_2_main;
        int mrna_total_edit_dist_2_main;
        int mrna_es_supp;
        double mrna_es_ratio;
    };
    vector<Orfhit> ORFs;

    cout << "ORF assembly loop (using DFS)\n";

    unordered_map<int, int> edited_support_DFS;
    for(size_t i = 0; i < DFS_results.size(); i++)
    {
        if(DFS_results[i].supp > 1)
        {
            Orfhit o;
            floParam flp;
            flp.M = DFS_results[i];
            flp.start_w_atg = Program_Options.start_with_ATG;
            flp.is_complete = Program_Options.draw_orf_incomplete;
            taORF orf_ = Find_Longest_ORF(flp);
            int total_supp, edited_supp;
            DFS_results[i].supp = Recalc_Support(flp.M, All, total_supp, edited_supp).size();

            edited_support_DFS[i] = edited_supp;

            o.start = orf_.tlstart;
            o.end = orf_.tlend;
            o.aa = TranslateORF(orf_.orf);
            o.nuc = orf_.orf;
            o.supp = DFS_results[i].supp;
            o.id = i;

            if(Program_Options.orf_length_min > 0 &&
                o.nuc.size() < Program_Options.orf_length_min) continue;

            //if(o.supp < 4) continue;
            ORFs.push_back(o);
        }
    }

    if(ORFs.size() > 0)
    {
        cout << "Sorting ORFs\n";
        sort(ORFs.begin(), ORFs.end(), [](Orfhit a, Orfhit b) { return a.aa.size() > b.aa.size();} );

        ofstream ofile_orf_aa((Program_Options.output_prefix+"_ORFs_aa.fasta").c_str());
        ofstream ofile_orf_nc((Program_Options.output_prefix+"_ORFs_nuc.fasta").c_str());
        ofstream ofile_orf_mrna((Program_Options.output_prefix+"_ORFs_mrna.fasta").c_str());
        DumpTless(Reference_Tless[0], (Program_Options.output_prefix+"_TlRef.txt").c_str());
        auto tlref = Reference_Tless[0];


        auto& main_orf = ORFs[Program_Options.main_orf_id];
        auto& main_orf_mp = DFS_results[main_orf.id];

        /*
            Calc revision stats
        */
        ofstream stat_file_es("Global_stats_AES_T-Aligner-rev.txt", std::ios_base::app);
        map<int, int> es_counter_good; int es_counter_total = 0;
        for(int i = 0; i < All.size(); i++)
        {
            int alt_states = Count_Alt_Edited_Sites(All[i], main_orf_mp, Reference_Tless[0]);
            es_counter_total++;
            for(int z = 1; z <= 20; z++)
            {
                if(alt_states >= z && alt_states <= 80)
                {
                    es_counter_good[z]++;
                }
            }
        }
        for(int z = 1; z <= 20; z++)
        {
            auto file_fa = Tokenize(Tokenize(Program_Options.fasta_genes_files[0], '/').back(), '.')[0];
            auto file_fq = Tokenize(Tokenize(Program_Options.fastq_input_files[0], '/').back(), '.')[0];
            stat_file_es << file_fa << "\t" << file_fq << "\t" << z << "\t" << es_counter_good[z] << "\t" << 100.0*es_counter_good[z]/es_counter_total << "\n";
            cout << file_fa << "\t" << file_fq << "\t" << z << "\t" << es_counter_good[z] << "\t" << 100.0*es_counter_good[z]/es_counter_total << "\n";
        }

        //return 0;
        vector<int> main_orf_edits(tlref.T.size(), -1000);
        int mo_begin = main_orf_mp.rf_start + main_orf.start;
        int mo_end = main_orf_mp.rf_start + main_orf.end;

        for(int i = main_orf.start; i < main_orf.end; i++)
        {
            main_orf_edits[main_orf_mp.rf_start+i] =
                main_orf_mp.mp.dT[i] - tlref.dT[main_orf_mp.rf_start+i];
        }

        for(int i = 0; i < ORFs.size(); i++)
        {
            auto& corf = ORFs[i];
            auto& corfmp = DFS_results[corf.id];

            ORFs[i].total_edit_dist_2_main = Get_Total_ES_dist(main_orf_mp, corfmp,
                Reference_Tless, 0, All, ORFs[i].es_ratio, ORFs[i].es_supp);

            ORFs[i].mrna_total_edit_dist_2_main = Get_Total_ES_dist_mRNA(main_orf_mp, corfmp,
                Reference_Tless, 0, All, ORFs[i].mrna_es_ratio, ORFs[i].mrna_es_supp);

            int edit_dist = 0;
            vector<int> current_orf_edits(tlref.T.size(), -1000);
            for(int j = corf.start; j < corf.end; j++)
            {
                current_orf_edits[corfmp.rf_start+j+1] =
                    corfmp.mp.dT[j+1] - tlref.dT[corfmp.rf_start+j+1];
                if(corfmp.rf_start+j < mo_begin || corfmp.rf_start+j >= mo_end)
                    current_orf_edits[corfmp.rf_start+j+1] = 0;
            }
            for(int j = 0; j < current_orf_edits.size(); j++)
            {
                if(current_orf_edits[j] == 0 ||
                    main_orf_edits[j]==-1000 || current_orf_edits[j]==-1000) continue;
                if(current_orf_edits[j] != main_orf_edits[j]) edit_dist++;
            }
            ORFs[i].edit_dist_2_main = edit_dist;

            // and same for mRNA
            vector<int> mrna_main_orf_edits(tlref.T.size(), -1000);

            for(int i = 1; i < main_orf_mp.mp.dT.size()-1; i++)
            {
                main_orf_edits[main_orf_mp.rf_start+i] =
                    main_orf_mp.mp.dT[i] - tlref.dT[main_orf_mp.rf_start+i];
            }
            int mrna_edit_dist = 0;
            vector<int> mrna_current_orf_edits(tlref.T.size(), -1000);
            for(int j = 1; j < corfmp.mp.dT.size()-1; j++)
            {
                mrna_current_orf_edits[corfmp.rf_start+j+1] =
                    corfmp.mp.dT[j+1] - tlref.dT[corfmp.rf_start+j+1];
                // no need to check anything out of mRNA's bounds
            }
            for(int j = 0; j < current_orf_edits.size(); j++)
            {
                if(current_orf_edits[j] == 0 ||
                    main_orf_edits[j]==-1000 || mrna_current_orf_edits[j]==-1000) continue;
                if(mrna_current_orf_edits[j] != main_orf_edits[j]) mrna_edit_dist++;
            }
            ORFs[i].mrna_edit_dist_2_main = mrna_edit_dist;
        }


        int coverage_frame_height = 0;
        vector<MappedPart> goodORFs;
        for(int i = 0; i < ORFs.size(); i++)
        {
            auto& that = DFS_results[ORFs[i].id];
            goodORFs.push_back(that);
            int mlen = PrintTless(that.mp).size();
            int Edits = 0; int Ins = 0; int Dels = 0;
            for(int j = ORFs[i].start; j <= ORFs[i].end; j++)
            {
                int Dindel = that.mp.dT[j+1] - tlref.dT[that.rf_start+j+1];
                if(Dindel != 0) Edits++;
                if(Dindel > 0) Ins++;
                if(Dindel < 0) Dels++;
            }
            for(int j = 0; j < Program_Options.draw_these_orfs.size(); j++)
            {
                if(Program_Options.draw_these_orfs[j] == i) coverage_frame_height += edited_support_DFS[ORFs[i].id];
            }

            ofile_orf_aa << ">ORF,ID=" << i << ",supp_total="<< ORFs[i].supp << ",supp_edited=" <<
                edited_support_DFS[ORFs[i].id] << ",mRNA_length=" <<
                mlen << ",Edits=" << Edits << ",(inserts:" <<
                Ins << ",dels:" << Dels << "),Dist=" <<
                ORFs[i].edit_dist_2_main << ",ES_Dist=" <<
                ORFs[i].total_edit_dist_2_main << ",ES_ratio=" << ORFs[i].es_ratio  << ",ES_Supp=" <<
                ORFs[i].es_supp << ",mRNA-Dist=" <<
                ORFs[i].mrna_edit_dist_2_main << ",mRNA-ES_Dist=" <<
                ORFs[i].mrna_total_edit_dist_2_main << ",mRNA-ES_ratio=" <<
                ORFs[i].mrna_es_ratio << ",mRNA-ES_Supp=" <<
                ORFs[i].mrna_es_supp << "\n" << ORFs[i].aa  << "\n";
            ofile_orf_nc << ">ORF,ID=" << i << ",supp_total="<< ORFs[i].supp << ",supp_edited=" <<
                edited_support_DFS[ORFs[i].id] << ",mRNA_length=" <<
                mlen << ",Edits=" << Edits << ",(inserts:" <<
                Ins << ",dels:" << Dels << "),Dist=" <<
                ORFs[i].edit_dist_2_main << ",ES_Dist=" <<
                ORFs[i].total_edit_dist_2_main << ",ES_ratio=" << ORFs[i].es_ratio << ",ES_Supp=" <<
                ORFs[i].es_supp << ",mRNA-Dist=" <<
                ORFs[i].mrna_edit_dist_2_main << ",mRNA-ES_Dist=" <<
                ORFs[i].mrna_total_edit_dist_2_main << ",mRNA-ES_ratio=" <<
                ORFs[i].mrna_es_ratio << ",mRNA-ES_Supp=" <<
                ORFs[i].mrna_es_supp << "\n" << ORFs[i].nuc  << "\n";

            ofile_orf_mrna << ">ORF,ID=" << i << ",supp_total="<< ORFs[i].supp << ",supp_edited=" <<
                edited_support_DFS[ORFs[i].id] << ",mRNA_length=" <<
                mlen << ",Edits=" << Edits << ",(inserts:" <<
                Ins << ",dels:" << Dels << "),Dist=" <<
                ORFs[i].edit_dist_2_main << ",ES_Dist=" <<
                ORFs[i].total_edit_dist_2_main << ",ES_ratio=" << ORFs[i].es_ratio << ",ES_Supp=" <<
                ORFs[i].es_supp << ",mRNA-Dist=" <<
                ORFs[i].mrna_edit_dist_2_main << ",mRNA-ES_Dist=" <<
                ORFs[i].mrna_total_edit_dist_2_main << ",mRNA-ES_ratio=" <<
                ORFs[i].mrna_es_ratio << ",mRNA-ES_Supp=" <<
                ORFs[i].mrna_es_supp << "\n" << PrintTless(that.mp)  << "\n";

        }
        /*
        * Here is the place to start drawing frames
        * because supp values now influence the frame
        * height
        */
        coverage_frame_height *= 2.5;

        //SFP.frame_heights.back() = coverage_frame_height + 60;
        Draw_Frames(SFP);

        // always draw TP
        if(Program_Options.mode_lib_compare)
        {
            ADO_cloud.mColor = QColor(140,140,140,25);
            Draw_Cloudmatrix(matrix_lib2, Reference_Tless, 0, ADO_cloud);
            /*
            auto ADO_COMP = ADO_cloud;
            ADO_COMP.Scene = Scene; ADO_COMP.UpperCorner += QPointF(0, 40.0);
            floParam flp;
            flp.M = DFS_results[main_orf.id];
            flp.start_w_atg = Program_Options.start_with_ATG;
            flp.is_complete = Program_Options.draw_orf_incomplete;
            Draw_ORF(flp, Reference_Tless, 0, ADO_COMP);
            */
        }
        else Draw_TranslatabilityPlot2(matrix, goodORFs, Reference_Tless, 0, ADO_cloud);
        /*
        if(Program_Options.draw_translatability_plot == true)
        {
            Draw_TranslatabilityPlot2(matrix, goodORFs, Reference_Tless, 0, ADO_cloud);
        }
        */

        ADO_orf.Scene = Scene; ADO_orf.UpperCorner = SFP.GetFrameCenter(5) + QPointF(0, 40.0);
        ADO_orf.mPenWidth = 1;
        ADO_orf_cloud = ADO_orf; ADO_orf_cloud.mColor = QColor(140,140,140,25);
        //ADO_orf_cloud.UpperCorner += QPointF(0, -10.0);
        ADO_orf_tp = ADO_orf; ADO_orf_tp.UpperCorner = SFP.GetFrameCenter(3) + QPointF(0, 40.0);
        auto ADO_ORF_supp = ADO_orf;

        ADO_ORF_supp.UpperCorner = upperleft_point + QPointF(0, 210);
        ADO_ORF_supp.Scene = MappedReadsScene;
        QGraphicsRectItem* whiteBackgound = new QGraphicsRectItem(
            0, 0, SFP.Reference.T.size() * SFP.x_axis_step + 100, coverage_frame_height+200);
        whiteBackgound->setBrush(QBrush(QColor(255,255,255,255)));
        QPen bbPen = QPen(QColor(0,0,0,255)); bbPen.setWidth(4);
        whiteBackgound->setPen(bbPen);
        ADO_ORF_supp.Scene->addItem(whiteBackgound);


        /*
        if(Program_Options.btlc_draw)
            ADO_ORF_supp.UpperCorner = SFP.GetFrameCenter(7) - QPointF(0, coverage_frame_height / 2);
        else
            ADO_ORF_supp.UpperCorner = SFP.GetFrameCenter(6) - QPointF(0, coverage_frame_height / 2);
        */

        Draw_Cloudmatrix(matrix, Reference_Tless, 0, ADO_orf_cloud);
        ADO_orf_cloud.Scene = MappedReadsScene;
        ADO_orf_cloud.UpperCorner = QPointF(0, 100) + QPointF(0, 40.0);
        Draw_Cloudmatrix(matrix, Reference_Tless, 0, ADO_orf_cloud);


        for(int i = 0; i < Program_Options.draw_these_orfs.size(); i++)
        {
            ADO_orf.mColor = QColor(255, 0, 0, 255); // default color is RED
            ADO_orf.UpperCorner += QPointF(0, i*(i%2==1?1:-1)*2);
            ADO_orf_tp.UpperCorner += QPointF(0, i*(i%2==1?1:-1)*2);
            ADO_orf_cloud.UpperCorner += QPointF(0, i*(i%2==1?1:-1)*2);

            if(i < Program_Options.draw_orfs_color_rgb.size())
            {
                ADO_orf.mColor = QColor(Program_Options.draw_orfs_color_rgb[i][0],
                    Program_Options.draw_orfs_color_rgb[i][1],
                    Program_Options.draw_orfs_color_rgb[i][2], 255);

            }
            ADO_orf_tp.mColor = ADO_orf.mColor;
            ADO_ORF_supp.mColor = ADO_orf.mColor;
            ADO_orf_cloud.mColor = ADO_orf.mColor;
            floParam flp;
            flp.M = DFS_results[ORFs[Program_Options.draw_these_orfs[i]].id];
            flp.start_w_atg = Program_Options.start_with_ATG;
            flp.is_complete = Program_Options.draw_orf_incomplete;
            Draw_ORF(flp, Reference_Tless, 0, ADO_orf);
            //if(Program_Options.draw_translatability_plot == true)
            //{
            Draw_ORF(flp, Reference_Tless, 0, ADO_orf_tp);
            //}
            Draw_ORF_Mapped_Reads(flp.M, All, Reference_Tless, 0, ADO_ORF_supp);
            Draw_ORF(flp, Reference_Tless, 0, ADO_orf_cloud);

            auto ADO_Aes = ADO_orf; ADO_Aes.UpperCorner = SFP.GetFrameCenter(4) + QPointF(0, 10*i -5.0 -5.0*Program_Options.draw_these_orfs.size());
            Draw_Alternatively_Edited_Sites(DFS_results[main_orf.id], flp.M, Reference_Tless, 0, ADO_Aes);
        }

        if(Program_Options.btlc_draw)
        {
            ADO_orf.mColor = QColor(255, 0, 0, 255); // default color is RED
            auto ADO_BTLC = ADO_orf;
            ADO_BTLC.Scene = Scene; ADO_BTLC.UpperCorner = SFP.GetFrameCenter(6) + QPointF(0, 40.0);
            Program_Options.btlc_pos -=2;
            floParam flp;
            flp.M = DFS_results[main_orf.id];
            flp.start_w_atg = Program_Options.start_with_ATG;
            flp.is_complete = Program_Options.draw_orf_incomplete;

            Draw_ORFs_BTLC(matrix, goodORFs, Reference_Tless, 0, ADO_BTLC, Program_Options);

            Draw_ORF(flp, Reference_Tless, 0, ADO_BTLC);
        }

        else if(Program_Options.mode_lib_compare)
        {
            ADO_orf.mColor = QColor(255, 0, 0, 255); // default color is RED
            auto ADO_COMP = ADO_orf;
            ADO_COMP.Scene = Scene; ADO_COMP.UpperCorner = SFP.GetFrameCenter(6) + QPointF(0, 40.0);
            floParam flp;
            flp.M = DFS_results[main_orf.id];
            flp.start_w_atg = Program_Options.start_with_ATG;
            flp.is_complete = Program_Options.draw_orf_incomplete;

            Draw_ORFs_Compare_Libs(matrix, matrix_lib2, flp, Reference_Tless, 0, ADO_COMP, Program_Options);

            Draw_ORF(flp, Reference_Tless, 0, ADO_COMP);
        }
    }
    else
    {
        cout << "No ORFs to draw, skipping the stage...\n";
        Draw_Frames(SFP);

    }

    Draw_Dotmatrix(matrix, Reference_Tless, 0, ADO_dots);
    map<int, vector<float> > CoverageHist;
    Draw_Coverage_Histogramm(matrix, Reference_Tless, 0, ADO_cover, CoverageHist);

    // dump coverage hist
    ofstream coverhist_out((Program_Options.output_prefix+"_histogram.txt").c_str());
    for(auto& x : CoverageHist) coverhist_out << x.first << "\t" <<
        x.second[0] << "\t" << x.second[1] << "\t" << x.second[2] << "\t" << "\n";

    // EARLY EXIT
    //return 0;

    cout << "\nRendering scene to PNG\nPlease wait...\n\n";
    QImage jpeg_image(20000,20000, QImage::Format_ARGB32);
    QPainter jpeg_painter(&jpeg_image);
    Scene->render(&jpeg_painter);
    jpeg_painter.end();
    jpeg_image.save((Program_Options.output_prefix+"_Image.png").c_str(), "PNG");
    //jpeg_image.save((Program_Options.output_prefix+"_Image.tiff").c_str(), "TIFF");
    //jpeg_image.save((Program_Options.output_prefix+"_Image.jpg").c_str(), "JPG");

    QImage jpeg_image_cover(20000,20000, QImage::Format_ARGB32);
    QPainter jpeg_painter_cover(&jpeg_image_cover);
    MappedReadsScene->render(&jpeg_painter_cover);
    jpeg_painter_cover.end();
    jpeg_image_cover.save((Program_Options.output_prefix+"_Image_cover.png").c_str(), "PNG");
    //jpeg_image_cover.save((Program_Options.output_prefix+"_Image_cover.tiff").c_str(), "TIFF");
    //jpeg_image_cover.save((Program_Options.output_prefix+"_Image_cover.jpg").c_str(), "JPG");


    QImage jpeg_image_low(5000,5000, QImage::Format_ARGB32);
    QPainter jpeg_painter_low(&jpeg_image_low);
    Scene->render(&jpeg_painter_low);
    jpeg_painter_low.end();
    jpeg_image_low.save((Program_Options.output_prefix+"_Image_lowres.png").c_str(), "PNG");
    //jpeg_image_low.save((Program_Options.output_prefix+"_Image_lowres.tiff").c_str(), "TIFF");
    //jpeg_image_low.save((Program_Options.output_prefix+"_Image_lowres.jpg").c_str(), "JPG");

    if(Program_Options.print_pdf)
    {
        cout << "\nRendering scene to PDF...\n";
        QPrinter printer;
        printer.setOutputFileName((Program_Options.output_prefix+"_Image.pdf").c_str());
        printer.setOutputFormat( QPrinter::PdfFormat );
        //printer.setResolution(100);
        printer.setPageSize(QPrinter::A4);
        printer.setOrientation(QPrinter::Landscape);
        QPainter painter(&printer);
        Scene->render(&painter);
        painter.end();
    }
    /*
    if(Program_Options.print_svg)
    {
        cout << "\nRendering scene to SVG...\n";
        QRectF newSceneRect;
        QGraphicsScene* tempScene = new QGraphicsScene(Scene->sceneRect());
        tempScene->setBackgroundBrush(QBrush(Qt::transparent));
        tempScene->setItemIndexMethod(QGraphicsScene::BspTreeIndex);
        foreach(QGraphicsItem* item, Scene->items()) {
            newSceneRect |= item->mapToScene(item->boundingRect()).boundingRect();
            tempScene->addItem(item);
        }
        tempScene->setSceneRect(newSceneRect);
        tempScene->clearSelection();
        QSize sceneSize = newSceneRect.size().toSize();

        QSvgGenerator generator;
        generator.setFileName((Program_Options.output_prefix+"_Image.svg").c_str());
        generator.setSize(sceneSize);
        generator.setViewBox(QRect(0, 0, sceneSize.width(), sceneSize.height()));
        QPainter painter;
        painter.begin(&generator);
        tempScene->render(&painter);
        painter.end();
        delete tempScene;
    }
    */
    cout << "All done! T-Aligner " << g_taligner_version << "\n";
    return 0;
}
