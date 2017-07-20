// by J@ :)

#ifndef FRAMESDRAWER_H
#define FRAMESDRAWER_H

#include "dna_util.h"
#include <QtGui>

using namespace std;
using namespace libdna;

struct SceneFramesParam
{
    QGraphicsScene* scene;
    QPointF UpperCorner;
    int n_frames;
    double x_axis_step;
    vector<int> frame_heights;
    vector<string> frame_names;
    TlessDNA Reference;
	string mapping_stats;
    string heading_label;
    QPixmap pixmap;

    QPointF GetFrameCenter(int frame)
    {
        double dy = UpperCorner.y()+80; int i = 0;
        for(; i < frame-1; i++)
        {
            dy += frame_heights[i];
            dy += 0;
        }
        dy += frame_heights[i] / 2;
        QPointF point(UpperCorner.x(), dy);
        return point;
    }
};


void Draw_Reference_Coordinates(SceneFramesParam& SFP, QPointF where)
{
    int frame_width =
        SFP.Reference.T.size() * SFP.x_axis_step + 100;

    string ref_sequence = PrintTless(SFP.Reference);

    QFont MarksFont("Courier New", 8, QFont::Bold);
    QPen MarksPen;
    MarksPen.setColor(QColor(0,0,0));
    MarksPen.setWidth(1);
    QPen InsPen;
    InsPen.setColor(QColor(250,0,0));
    InsPen.setWidth(1);

    int real_seq_pos = SFP.Reference.dT[0];
    for(int dx = 0; dx < SFP.Reference.T.size(); dx++)
    {
        if(dx % 10 == 0 && dx >= 10)
        {
            auto MarkText = SFP.scene->addText(
                i2s(real_seq_pos).c_str(), MarksFont);
            MarkText->setPos(where.x()+(dx)*SFP.x_axis_step, where.y()+10);
        }
        real_seq_pos++;
        if(dx < SFP.Reference.dT.size())
        {
            if(dx >= 10)
            for(int i = 0; i < SFP.Reference.dT[dx]; i++)
            {
                QGraphicsRectItem* Tinsert = new QGraphicsRectItem(
                    where.x()+(dx)*SFP.x_axis_step, where.y()+25+2*i, SFP.x_axis_step, 1);
                Tinsert->setPen(InsPen);
                SFP.scene->addItem(Tinsert);
            }
            real_seq_pos += SFP.Reference.dT[dx];
        }
    }
}

void Draw_Frames(SceneFramesParam& SFP)

{
    // Black bold pen for frames
    QPen FramePen;
    FramePen.setColor(QColor(0,0,0));
    FramePen.setWidth(4);

    int frame_width =
        SFP.Reference.T.size() * SFP.x_axis_step + 100;

    // Black Times New Roman 12 font for labels
    QFont InfoFont("Times New", 12, QFont::Bold);
    QFont InfoFontL("Times New", 20, QFont::Bold);

    QGraphicsRectItem* StatsFrame = new QGraphicsRectItem(
        SFP.UpperCorner.x(), SFP.UpperCorner.y(), frame_width, 50);
    //StatsFrame->setBrush(QBrush(QColor(240,230,255,100)));
    // make white background
    StatsFrame->setBrush(QBrush(QColor(255,255,255,255)));
    StatsFrame->setPen(FramePen);
    SFP.scene->addItem(StatsFrame);

    auto HeadingText = SFP.scene->addText(
        SFP.heading_label.c_str(), InfoFont);
    HeadingText->setPos(SFP.UpperCorner.x()+10, SFP.UpperCorner.y()+5);

    auto MappingText = SFP.scene->addText(
        SFP.mapping_stats.c_str(), InfoFont);
    MappingText->setPos(SFP.UpperCorner.x()+10, SFP.UpperCorner.y()+20);

    QPointF coord_point(SFP.UpperCorner.x(), SFP.UpperCorner.y()+40);

    // 100

    int current_frame_Y = SFP.UpperCorner.y()+80;

    // make white background
    QGraphicsRectItem* whiteBackgound = new QGraphicsRectItem(
        SFP.UpperCorner.x(), SFP.UpperCorner.y() + 40, frame_width, 40);
    whiteBackgound->setBrush(QBrush(QColor(255,255,255,255)));
    whiteBackgound->setPen(FramePen);
    SFP.scene->addItem(whiteBackgound);
    Draw_Reference_Coordinates(SFP, coord_point);

    for(size_t frame = 0; frame < SFP.n_frames; frame++)
    {
        QGraphicsRectItem* Frame = new QGraphicsRectItem(
            SFP.UpperCorner.x(), current_frame_Y, frame_width,
            SFP.frame_heights[frame]);
        Frame->setBrush(QBrush(QColor(255,255,255,255)));
        Frame->setPen(FramePen);
        SFP.scene->addItem(Frame);

        QGraphicsPixmapItem* logoItem = new QGraphicsPixmapItem(SFP.pixmap);
        logoItem->scale(0.1,0.1);
        logoItem->setPos(frame_width-50,current_frame_Y+10);
        SFP.scene->addItem(logoItem);

        auto FrameHeadingText = SFP.scene->addText(
            SFP.frame_names[frame].c_str(), InfoFontL);
        FrameHeadingText->setPos(SFP.UpperCorner.x() + frame_width - 210, current_frame_Y+SFP.frame_heights[frame] - 50);
        FrameHeadingText->setZValue(100);
        current_frame_Y += SFP.frame_heights[frame];
        current_frame_Y += 0;
    }
}

#endif
