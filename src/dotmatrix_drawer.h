// by J@ :)

#ifndef DOTDRAWER_H
#define DOTDRAWER_H

#include <random>
#include <algorithm>
#include <utility>
#include <math.h>

#include "fastq_reader.h"
#include "orf_finder.h"
// For plotter
#include <QtGui>

using namespace std;
using namespace libdna;

static std::random_device rand_dev;
static std::mt19937 rand_gen(rand_dev());

const int max_indel_l = 20;
const int draw_dot_rdsupp_min = 2;
const double x_axis_step = 7.0;

struct AbstractDrawOptions
{
	QPointF UpperCorner;
	QGraphicsScene* Scene;
	QColor mColor;
	vector<QColor> vColor;
	QPen mPen;
	QBrush mBrush;
	int mPenWidth;
};

void Draw_Dotmatrix(vector<AlignedRow>& matrix,
	vector<TlessDNA>& Reference_base,
	size_t Reference_id,
	AbstractDrawOptions ADO)
{
    auto& rtless = Reference_base[Reference_id];

	map<int, vector<int> > sorted_alignments;

	for(size_t i = 0; i < matrix.size(); i++)
	{
		// Should plot alignments only for given reference
		if(matrix[i].ref == Reference_id)
		{
			sorted_alignments[matrix[i].rf_start].push_back(i);
		}
	}


	map<int, map<int,int> > Editing_Matrix;
	for(size_t i = 0; i < rtless.T.size(); i++)
	{
		if(sorted_alignments.count(i))
		{
			for(auto& id : sorted_alignments[i])
			{
				auto& alignment = matrix[id];
				for(size_t p = 0; p < alignment.r.size(); p++)
					Editing_Matrix[alignment.r[p]][i+p]++;
			}
		}
	}


	map<int, map<int,int> > Actual_dot_alpha;
	for(int i = 0; i < rtless.T.size(); i++)
    {
		int current_col_summ = 0; int current_col_max = 0;
		for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
			current_col_summ += Editing_Matrix[k][i];
			if(Editing_Matrix[k][i] > current_col_max)
				current_col_max = Editing_Matrix[k][i];
		}
		for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
			Actual_dot_alpha[k][i] = 500*(1.0*Editing_Matrix[k][i] / current_col_max);
		}
	}

    // Draw matrix cycle
    double pos_x = 0;
    for(int i = 0; i < rtless.T.size(); i++)
    {
        pos_x += x_axis_step;
        for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
            int cell_rd_count = Editing_Matrix[k][i];
            if(cell_rd_count >= draw_dot_rdsupp_min)
            {
                int alpha = min(255.0, 25.5+Actual_dot_alpha[k][i]);
				int red   = 0;
				int green = 0;
				int blue  = max(0, Actual_dot_alpha[k][i]-250);

				QGraphicsEllipseItem* MatrixDot = new QGraphicsEllipseItem(
					pos_x + 1.0, ADO.UpperCorner.y() + x_axis_step*(-1)*k + 1.0, 6.0, 6.0);
				MatrixDot->setPen(QPen(QColor(red, green, blue, alpha)));
				MatrixDot->setBrush(QBrush(QColor(red, green, blue, alpha)));
				ADO.Scene->addItem(MatrixDot);

            }
        }
    }

    // legend
    QGraphicsRectItem* legendBar = new QGraphicsRectItem(
		100.0, ADO.UpperCorner.y() + 80, x_axis_step*2, x_axis_step*1.5);
    legendBar->setBrush(QBrush(QColor(255,255,255,255)));
    ADO.Scene->addItem(legendBar);
    for(double i = 0.0; i < 100.0; i+= 5.0)
    {
        QGraphicsRectItem* legendBar = new QGraphicsRectItem(
			100.0+2*x_axis_step*i/5.0+2*x_axis_step,
            ADO.UpperCorner.y() + 80,
            x_axis_step*2,
            x_axis_step*1.5);
        legendBar->setBrush(QBrush(QColor(0,0,max(0.0, 500*(i/100)-250), min(255.0, 25.5+500*(i/100)))));
        ADO.Scene->addItem(legendBar);
    }
    QFont MarksFont("Courier New", 15, QFont::Bold);
    auto legMin = ADO.Scene->addText("0", MarksFont);
    legMin->setPos(85.0, ADO.UpperCorner.y() + 80);
    auto legMax = ADO.Scene->addText("100%", MarksFont);
    legMax->setPos(100.0+42*x_axis_step+5.0, ADO.UpperCorner.y() + 80);

}

void Draw_Coverage_Histogramm(vector<AlignedRow>& matrix,
	vector<TlessDNA>& Reference_base,
	size_t Reference_id,
	AbstractDrawOptions ADO,
	map<int, vector<float> >& histogram)
{

	map<int, map<string, double> > editing_data;
	int max_col_value = 0;
    for(auto& alignment : matrix)
	{
		for(int i = 0; i < alignment.r.size(); i++)
		{
			if(alignment.r[i] == 0) editing_data[i+alignment.rf_start]["ref"]++;
			else if(alignment.r[i] > 0) editing_data[i+alignment.rf_start]["ins"]++;
			else editing_data[i+alignment.rf_start]["del"]++;
		}
	}

	// determine maximum val
	for(int rp = 0; rp < Reference_base[Reference_id].T.size(); rp++)
	{
		if(editing_data.count(rp)==0) {
			editing_data[rp]["ref"] = 0;
			editing_data[rp]["ins"] = 0;
			editing_data[rp]["del"] = 0;
		}
		else {
			if(editing_data[rp].count("ref") == 0) editing_data[rp]["ref"] = 0;
			if(editing_data[rp].count("ins") == 0) editing_data[rp]["ins"] = 0;
			if(editing_data[rp].count("del") == 0) editing_data[rp]["del"] = 0;
		}
		int cur_col_value = editing_data[rp]["ref"] +
							editing_data[rp]["ins"] +
							editing_data[rp]["del"];
		if(cur_col_value > max_col_value) max_col_value = cur_col_value;
	}


	// draw
	double by = ADO.UpperCorner.y()-90.0;
	const double h_size_max = 180.0;
	QBrush refBrush = QBrush(QColor("#41D167"));
	QBrush insBrush = QBrush(QColor("#4163D1"));
	QBrush delBrush = QBrush(QColor("#D141AB"));


    QPen EditingGraphPen;
	EditingGraphPen.setWidth(1.0);
	EditingGraphPen.setColor(QColor(30, 30, 30, 255));
    QPen EditingGraphPenSolid = EditingGraphPen;
    EditingGraphPen.setStyle(Qt::DotLine);

	QPen nulPen = QPen(Qt::color0);
	for(int rp = 0; rp < Reference_base[Reference_id].T.size(); rp++)
	{
		double c_bar_y  = 0.0;
		double d_y = max(1.0, editing_data[rp]["ref"]*h_size_max/max_col_value);
		if(editing_data[rp]["ref"]==0) d_y=0.0;
		QGraphicsRectItem* refBar = new QGraphicsRectItem(
			x_axis_step*(1.0+rp), by +c_bar_y, x_axis_step, d_y);
		c_bar_y += d_y;
		d_y = max(1.0, editing_data[rp]["ins"]*h_size_max/max_col_value);
		if(editing_data[rp]["ins"]==0) d_y=0.0;
		QGraphicsRectItem* insBar = new QGraphicsRectItem(
			x_axis_step*(1.0+rp), by +c_bar_y, x_axis_step, d_y);
		c_bar_y += d_y;
		d_y = max(1.0, editing_data[rp]["del"]*h_size_max/max_col_value);
		if(editing_data[rp]["del"]==0) d_y=0.0;
		QGraphicsRectItem* delBar = new QGraphicsRectItem(
			x_axis_step*(1.0+rp), by +c_bar_y, x_axis_step, d_y);
		refBar->setPen(nulPen); refBar->setBrush(refBrush);
		insBar->setPen(nulPen); insBar->setBrush(insBrush);
		delBar->setPen(nulPen); delBar->setBrush(delBrush);

		ADO.Scene->addItem(refBar);
		ADO.Scene->addItem(insBar);
		ADO.Scene->addItem(delBar);
		int summ = editing_data[rp]["ref"] + editing_data[rp]["ins"] + editing_data[rp]["del"];
		histogram[rp].push_back(editing_data[rp]["ref"]/summ);
		histogram[rp].push_back(editing_data[rp]["ins"]/summ);
		histogram[rp].push_back(editing_data[rp]["del"]/summ);

        double editing_pp = by + (h_size_max*(editing_data[rp]["ins"]+editing_data[rp]["del"]) / summ);
        if(editing_pp - by > 0.05 * h_size_max) {
            QGraphicsLineItem* editingPercent = new QGraphicsLineItem(x_axis_step*(1.0+rp),
                editing_pp, x_axis_step*(1.0+rp)+x_axis_step, editing_pp);
            editingPercent->setPen(EditingGraphPenSolid);
            ADO.Scene->addItem(editingPercent);
            QGraphicsLineItem* editingPercentBase = new QGraphicsLineItem(
                x_axis_step*(1.0+rp)+x_axis_step/2, by+2,
                x_axis_step*(1.0+rp)+x_axis_step/2, editing_pp);
            editingPercentBase->setPen(EditingGraphPen);
            ADO.Scene->addItem(editingPercentBase);
        }
	}


    // add Y axis
    QFont MarksFont("Courier New", 15, QFont::Bold);
    QFont MarksFontS("Courier New", 12, QFont::Bold);

    QGraphicsRectItem* OY = new QGraphicsRectItem(30.0, by, 2, h_size_max);
    OY->setPen(nulPen); OY->setBrush(QBrush(QColor(0,0,0,255)));
    ADO.Scene->addItem(OY);
    //QGraphicsRectItem* OYp = new QGraphicsRectItem(31.5, by, 1, h_size_max);
    //OYp->setPen(nulPen); OYp->setBrush(QBrush(QColor(255,0,0,255)));
    //ADO.Scene->addItem(OYp);

    for(int i = 0; i < 5; i++) {
        QGraphicsRectItem* Ytick = new QGraphicsRectItem(30.0, by + i * h_size_max / 4, 7.0, 2.0);
        Ytick->setPen(nulPen); Ytick->setBrush(QBrush(QColor(0,0,0,255)));
        ADO.Scene->addItem(Ytick);
        if(i < 4) {
            auto Yticktext = ADO.Scene->addText(QString(i2s(0.0 + i*max_col_value/4).c_str()) +
                " / " + QString(i2s(0 + i*100/4).c_str()) + "%", MarksFont);
            Yticktext->setPos(35.0, by - 5.0 + i * h_size_max / 4);
        }
    }
    auto Ytick100p = ADO.Scene->addText(QString(i2s(0.0+max_col_value).c_str()) + " / 100%", MarksFont);
        Ytick100p->setDefaultTextColor(QColor(0,0,0,255));
        Ytick100p->setPos(35.0, by - 12.0 + h_size_max);
    QGraphicsRectItem* legRef = new QGraphicsRectItem(
			200.0, by - 12.0 + h_size_max, 2*x_axis_step, 2*x_axis_step);
    QGraphicsRectItem* legIns = new QGraphicsRectItem(
			280, by - 12.0 + h_size_max, 2*x_axis_step, 2*x_axis_step);
    QGraphicsRectItem* legDel = new QGraphicsRectItem(
			360, by - 12.0 + h_size_max, 2*x_axis_step, 2*x_axis_step);

    auto legRefText = ADO.Scene->addText("Not\nedited", MarksFontS);
    legRefText->setPos(220, by - 21.0 + h_size_max);
    auto legInsText = ADO.Scene->addText("INS", MarksFont);
    legInsText->setPos(300, by - 20.0 + h_size_max);
    auto legDelText = ADO.Scene->addText("DEL", MarksFont);
    legDelText->setPos(380, by - 20.0 + h_size_max);


    legRef->setPen(nulPen); legRef->setBrush(refBrush);
	legIns->setPen(nulPen); legIns->setBrush(insBrush);
	legDel->setPen(nulPen); legDel->setBrush(delBrush);

	ADO.Scene->addItem(legRef);
	ADO.Scene->addItem(legIns);
	ADO.Scene->addItem(legDel);

}

void Draw_ORF(floParam& ORF,
	vector<TlessDNA>& Reference_base,
	size_t Reference_id,
	AbstractDrawOptions ADO
)
{

	QPen ORF_pen;
	ORF_pen.setWidth(ADO.mPenWidth);
	ORF_pen.setColor(ADO.mColor);

	QPainterPath ORF_path;
	auto O = Find_Longest_ORF(ORF);

	vector<int> R;
	for(int i = O.tlstart; i < O.tlend; i++) R.push_back(ORF.M.mp.dT[i+1]);

	for(size_t i = 0; i < R.size(); i++)
	{
		R[i] = R[i] -
			Reference_base[Reference_id].dT[ORF.M.rf_start+O.tlstart+1+i];
	}


	QPen MarkerPen = QPen(QColor(255,0,0));
	MarkerPen.setWidth(ADO.mPenWidth);
	QPolygonF StartCodonPolygon;

	QPointF BasePoint(0.5*x_axis_step+x_axis_step*(ORF.M.rf_start+O.tlstart),
		ADO.UpperCorner.y() + x_axis_step*(-1)*(ORF.M.mp.dT[O.tlstart] - Reference_base[Reference_id].dT[ORF.M.rf_start+O.tlstart]+1.0));

	//QPointF BasePoint(0.5*x_axis_step+x_axis_step*(ORF.M.rf_start+O.tlstart+1),
	//	ADO.UpperCorner.y() /*+4.0*x_axis_step*/ + x_axis_step*(-1)*(R[0]+1.0));
	StartCodonPolygon << BasePoint + QPointF(3,0) << BasePoint + QPointF(-1,-3) <<
		BasePoint + QPointF(-1,3);

	QGraphicsPolygonItem* ORFstart = new QGraphicsPolygonItem(StartCodonPolygon);
	ORFstart->setPen(MarkerPen);

	ORF_path.moveTo(BasePoint.x(), BasePoint.y());
	int u = 0;
	for(; u < R.size(); u++)
	{
		if(R[u] > 17 || R[u] < -17) continue;
		ORF_path.lineTo(x_axis_step*(u+0.5+ORF.M.rf_start+O.tlstart+1),
			ADO.UpperCorner.y() /*+4.0*x_axis_step*/ + x_axis_step*(-1)*(R[u]+1.0));
	}
	u--;
	if(R[u] > 17 || R[u] < -17) u--;
	QGraphicsRectItem* ORFend = new QGraphicsRectItem(x_axis_step*(u+0.5+ORF.M.rf_start+O.tlstart+1)-3,
		ADO.UpperCorner.y() /*+4.0*x_axis_step*/ + x_axis_step*(-1)*(R[u]+1.0)-3,6,6);
	ORFend->setPen(MarkerPen);

	QGraphicsPathItem* ORF_Item = new QGraphicsPathItem(ORF_path);
	ORF_Item->setPen(ORF_pen);
	ADO.Scene->addItem(ORF_Item);
	ADO.Scene->addItem(ORFstart);
	ADO.Scene->addItem(ORFend);
}


void Draw_Cloudmatrix(vector<AlignedRow>& matrix,
	vector<TlessDNA>& Reference_base,
	size_t Reference_id,
	AbstractDrawOptions ADO
)
{

	// draw axis

	for(int i = -17; i <= 7; i++)
	{
		QGraphicsLineItem* horizLine = new QGraphicsLineItem(
		x_axis_step,
		ADO.UpperCorner.y()+ 8.0 + x_axis_step*i,
		ADO.Scene->width()-2*x_axis_step,
		ADO.UpperCorner.y()+ 8.0 + x_axis_step*i);
		horizLine->setPen(QPen(QColor(0,0,0,200)));
		ADO.Scene->addItem(horizLine);
	}

	// draw cloud

    auto& rtless = Reference_base[Reference_id];
	uniform_real_distribution<> R(0.0, 0.01+x_axis_step/3);
	uniform_real_distribution<> Psi(-3.1415926, 3.1415926);
	map<double, map<double, int> > Data_cloud;

	for(auto& align : matrix)
	{
		for(size_t i = 0; i < align.r.size(); i++)
		{

			double tX = 1.0*x_axis_step*static_cast<double>(align.rf_start+i+1.5) +
			R(rand_gen)*qCos(Psi(rand_gen));

			double tY = ADO.UpperCorner.y() + //4.0*x_axis_step +
			1.0*x_axis_step*(-1)*static_cast<double>(align.r[i]+1) +
			R(rand_gen)*qSin(Psi(rand_gen));

			tX = floor( tX * 10 ) / 10;
			tY = floor( tY * 10 ) / 10;
			Data_cloud[tX][tY]++;

		}
	}
	int nomA = ADO.mColor.alpha();
	for(auto& pX : Data_cloud)
	{
		 for(auto& pY : pX.second)
		 {
			 double tX = pX.first;
			 double tY = pY.first;
			 int tA = pY.second;

			 QColor tC = ADO.mColor;
			 int cA = nomA*tA; if(cA > 255) cA = 255;
			 tC.setAlpha(cA);
			 QGraphicsEllipseItem* MatrixDot = new QGraphicsEllipseItem(
				 tX-2.0,
				 tY-2.0,
				 4.0, 4.0);
			 MatrixDot->setPen(QPen(tC));
			 MatrixDot->setBrush(QBrush(tC));
			 ADO.Scene->addItem(MatrixDot);
		 }
	}
}


void Draw_Cloudmatrix_v2(vector<AlignedRow>& matrix,
	vector<TlessDNA>& Reference_base,
	size_t Reference_id,
	AbstractDrawOptions ADO
)
{
    auto& rtless = Reference_base[Reference_id];
	// draw axis

	for(int i = -17; i <= 7; i++)
	{
		QGraphicsLineItem* horizLine = new QGraphicsLineItem(
		x_axis_step,
		ADO.UpperCorner.y()+ 8.0 + x_axis_step*i,
		ADO.Scene->width()-2*x_axis_step,
		ADO.UpperCorner.y()+ 8.0 + x_axis_step*i);
		horizLine->setPen(QPen(QColor(0,0,0,200)));
		ADO.Scene->addItem(horizLine);
	}


    map<int, vector<int> > sorted_alignments;
	for(size_t i = 0; i < matrix.size(); i++) {
		if(matrix[i].ref == Reference_id) {
			sorted_alignments[matrix[i].rf_start].push_back(i);
		}
	}

	map<int, map<int,int> > Editing_Matrix;
	for(size_t i = 0; i < rtless.T.size(); i++)	{
		if(sorted_alignments.count(i)) {
			for(auto& id : sorted_alignments[i]) {
				auto& alignment = matrix[id];
				for(size_t p = 0; p < alignment.r.size(); p++)
					Editing_Matrix[alignment.r[p]][i+p]++;
			}
		}
	}


	map<int, map<int,int> > Actual_dot_alpha;
    int global_max = 0;
	for(int i = 0; i < rtless.T.size(); i++)
    {
		int current_col_summ = 0; int current_col_max = 0;
		for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
			current_col_summ += Editing_Matrix[k][i];
			if(Editing_Matrix[k][i] > current_col_max)
				current_col_max = Editing_Matrix[k][i];
		}
        if(current_col_max > global_max) global_max = current_col_max;

	}

    for(int i = 0; i < rtless.T.size(); i++)
    {
        for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
			Actual_dot_alpha[k][i] = 500*(1.0*Editing_Matrix[k][i] / global_max);
		}
    }



    // Draw matrix cycle
    double pos_x = 0;
    for(int i = 0; i < rtless.T.size(); i++)
    {
        pos_x += x_axis_step;
        for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
            int cell_rd_count = Editing_Matrix[k][i];
            if(cell_rd_count >= draw_dot_rdsupp_min)
            {
                int alpha = min(255.0, 25.5+Actual_dot_alpha[k][i]);
				int red   = 0;
				int green = 0;
				int blue  = max(0, Actual_dot_alpha[k][i]-250);

				QGraphicsEllipseItem* MatrixDot = new QGraphicsEllipseItem(
					pos_x + 1.0, ADO.UpperCorner.y() + x_axis_step*(-1)*k + 1.0, 6.0, 6.0);
				MatrixDot->setPen(QPen(QColor(red, green, blue, alpha)));
				MatrixDot->setBrush(QBrush(QColor(red, green, blue, alpha)));
				ADO.Scene->addItem(MatrixDot);

            }
        }
    }
}


void Draw_TranslatabilityPlot2(vector<AlignedRow>& matrix,
	vector<MappedPart>& translations,
	vector<TlessDNA>& Reference_base,
	size_t Reference_id,
	AbstractDrawOptions ADO)
{
    auto& rtless = Reference_base[Reference_id];

	map<int, vector<int> > sorted_alignments;

	for(size_t i = 0; i < matrix.size(); i++)
	{
		// Should plot alignments only for given reference
		if(matrix[i].ref == Reference_id)
		{
			sorted_alignments[matrix[i].rf_start].push_back(i);
		}
	}

	map<int, vector<vector<int> > > sorted_translations;
    /*
	for(size_t i = 0; i < translations.size(); i++)
	{
		auto& tr = translations[i];
		vector<int> r;
		for(int j = 0; j < (tr.mp.dT.size()-1) &&
		tr.rf_start+j+1 < rtless.dT.size(); j++)
		{
			r.push_back(tr.mp.dT[j+1] - rtless.dT[tr.rf_start+j+1]);
		}
		if(r.size() > 5)
			sorted_translations[tr.rf_start].push_back(r);
	}
	*/
    for(size_t i = 0; i < translations.size(); i++)
	{
        floParam flp;
        flp.M = translations[i];
        flp.start_w_atg = false;
        flp.is_complete = false;
		auto O = Find_Longest_ORF(flp);
        vector<int> r;

        for(int k = O.tlstart; k <= O.tlend; k++) r.push_back(flp.M.mp.dT[k+1]);
		for(size_t k = 0; k < r.size(); k++) r[k] = r[k] - rtless.dT[flp.M.rf_start+O.tlstart+1+k];
        if(r.size() > 5)
			sorted_translations[flp.M.rf_start+O.tlstart].push_back(r);
    }

	map<int, map<int,int> > Translated_Matrix;
	map<int, map<int,int> > Editing_Matrix;
	for(size_t i = 0; i < rtless.T.size(); i++)
	{
		if(sorted_alignments.count(i))
		{
			for(auto& id : sorted_alignments[i])
			{
				auto& alignment = matrix[id];
				for(size_t p = 0; p < alignment.r.size(); p++)
					Editing_Matrix[alignment.r[p]][i+p]++;
			}
		}
		if(sorted_translations.count(i))
		{
			for(auto& tr : sorted_translations[i])
			{
				for(size_t p = 0; p < tr.size(); p++)
				{
					Translated_Matrix[tr[p]][i+p]++;
				}
			}
		}
	}

	map<int, map<int,int> > Actual_dot_alpha;
	for(int i = 0; i < rtless.T.size(); i++)
    {
		int current_col_summ = 0; int current_col_max = 0;
		for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
			current_col_summ += Editing_Matrix[k][i];
			if(Editing_Matrix[k][i] > current_col_max)
				current_col_max = Editing_Matrix[k][i];
		}
		for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
			Actual_dot_alpha[k][i] = 500*(1.0*Editing_Matrix[k][i] / current_col_max);
		}
	}

    // Draw matrix cycle
    double pos_x = 0;
    for(int i = 0; i < rtless.T.size(); i++)
    {
        pos_x += x_axis_step;
        for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
            int cell_rd_count = Editing_Matrix[k][i];
			int alpha = min(255.0, 25.5+Actual_dot_alpha[k][i]);
			int red   = 0;
			int green = 0;
			int blue  = max(0, Actual_dot_alpha[k][i]-250);

			if(Translated_Matrix[k][i] > 0)
			{
				/*
                QGraphicsEllipseItem* MatrixDotE = new QGraphicsEllipseItem(
					pos_x+3.0,
					ADO.UpperCorner.y() + x_axis_step*(-1)*k - 1.5*x_axis_step, // - 3.0,
					8.0, 8.0);
				*/
                QGraphicsRectItem* MatrixDotE = new QGraphicsRectItem(
					pos_x-1.0,
					ADO.UpperCorner.y() + x_axis_step*(-1)*k - 1.5*x_axis_step-1.0, // - 3.0,
					9.0, 9.0);
				QPen ePen = QPen(QColor(0, 250, 20, 125));
				ePen.setWidth(1);
				MatrixDotE->setPen(ePen);
				MatrixDotE->setBrush(QBrush(QColor(0, 250, 20, 125)));
				ADO.Scene->addItem(MatrixDotE);
			}

            if(cell_rd_count >= draw_dot_rdsupp_min)
            {
                QGraphicsEllipseItem* MatrixDotBackground = new QGraphicsEllipseItem(
					pos_x + 1.0,
					ADO.UpperCorner.y() + x_axis_step*(-1)*k - 1.5*x_axis_step + 1.0, // - 2.0,
					6.0, 6.0);
				MatrixDotBackground->setPen(QPen(QColor(255,255,255,255)));
				MatrixDotBackground->setBrush(QBrush(QColor(255,255,255,255)));
				ADO.Scene->addItem(MatrixDotBackground);

				QGraphicsEllipseItem* MatrixDot = new QGraphicsEllipseItem(
					pos_x + 1.0,
					ADO.UpperCorner.y() + x_axis_step*(-1)*k - 1.5*x_axis_step + 1.0, // - 2.0,
					6.0, 6.0);
				MatrixDot->setPen(QPen(QColor(red, green, blue, alpha)));
				MatrixDot->setBrush(QBrush(QColor(red, green, blue, alpha)));
				ADO.Scene->addItem(MatrixDot);
            }
        }
    }
}



void Draw_ORFs_BTLC(vector<AlignedRow>& matrix,
	vector<MappedPart>& translations,
	vector<TlessDNA>& Reference_base,
	size_t Reference_id,
	AbstractDrawOptions ADO,
    TAlignerOptions Program_Options)
{
    int btlc_pos = Program_Options.btlc_pos;
    int edited_k = Program_Options.btlc_k;
    int count = Program_Options.btlc_count;
    auto& rtless = Reference_base[Reference_id];

	map<int, vector<int> > sorted_alignments;

	for(size_t i = 0; i < matrix.size(); i++)
	{
		// Should plot alignments only for given reference
		if(matrix[i].ref == Reference_id)
		{
			sorted_alignments[matrix[i].rf_start].push_back(i);
		}
	}

	vector<vector<vector<int> > > sorted_translations;
	for(int i = 0; i < translations.size(); i++)
	{
		floParam flp;
        flp.M = translations[i];
        flp.start_w_atg = false;
        flp.is_complete = true;
		auto O = Find_Longest_ORF(flp);
		vector<vector<int> > R(4);
		R[0].push_back(flp.M.rf_start + O.tlstart);
		R[2].push_back(i);
		R[3].push_back(O.orf.size());
		for(int k = O.tlstart; k <= O.tlend; k++) R[1].push_back(flp.M.mp.dT[k+1]);
		for(size_t k = 0; k < R[1].size(); k++) R[1][k] = R[1][k] - rtless.dT[flp.M.rf_start+O.tlstart+1+k];
		sorted_translations.push_back(R);
	}


	sort(sorted_translations.begin(), sorted_translations.end(), [](vector<vector<int> > A, vector<vector<int> > B)
		{ return A[3][0] > B[3][0]; });

	vector<int> accepted_orfs;

	for(size_t i = 0; i < sorted_translations.size(); i++)
	{
		auto& R = sorted_translations[i];
		for(int j = 0; j <  R[1].size(); j++)
		{
			int state = R[1][j];
			int tlp = R[0][0]+j;
			if(state == edited_k && tlp == btlc_pos) {
				accepted_orfs.push_back(i);
				break;
			}
		}
	}

	cout << "Computing ORFs for point (" << btlc_pos << ", " << edited_k <<  ")\n";

	map<int, map<int,int> > Editing_Matrix;
	for(size_t i = 0; i < rtless.T.size(); i++)
	{
		if(sorted_alignments.count(i)) {
			for(auto& id : sorted_alignments[i]) {
				auto& alignment = matrix[id];
				for(size_t p = 0; p < alignment.r.size(); p++)
					Editing_Matrix[alignment.r[p]][i+p]++;
			}
		}
	}

	map<int, map<int,int> > Actual_dot_alpha;
	for(int i = 0; i < rtless.T.size(); i++)
    {
		int current_col_summ = 0; int current_col_max = 0;
		for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
			current_col_summ += Editing_Matrix[k][i];
			if(Editing_Matrix[k][i] > current_col_max)
				current_col_max = Editing_Matrix[k][i];
		}
		for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
			Actual_dot_alpha[k][i] = 500*(1.0*Editing_Matrix[k][i] / current_col_max);
		}
	}

    // Draw matrix cycle
    double pos_x = 0;
    for(int i = 0; i < rtless.T.size(); i++)
    {
        pos_x += x_axis_step;
        for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
            int cell_rd_count = Editing_Matrix[k][i];
			int alpha = min(255.0, 25.5+Actual_dot_alpha[k][i]);
			int red   = 0;
			int green = 0;
			int blue  = max(0, Actual_dot_alpha[k][i]-250);

            if(cell_rd_count >= draw_dot_rdsupp_min)
            {
				QGraphicsEllipseItem* MatrixDot = new QGraphicsEllipseItem(
					pos_x + 1.0,
					ADO.UpperCorner.y() + x_axis_step*(-1)*k - 1.5*x_axis_step + 1.0, // - 2.0,
					6.0, 6.0);
				MatrixDot->setPen(QPen(QColor(red, green, blue, alpha)));
				MatrixDot->setBrush(QBrush(QColor(red, green, blue, alpha)));
				ADO.Scene->addItem(MatrixDot);
            }

			if(k == edited_k && i == btlc_pos)
			{
				QGraphicsEllipseItem* MatrixDot = new QGraphicsEllipseItem(
					pos_x + 1.0,
					ADO.UpperCorner.y() + x_axis_step*(-1)*k - 1.5*x_axis_step + 1.0, // - 2.0,
					6.0, 6.0);
				MatrixDot->setPen(QPen(QColor(255,0,0,255)));
				MatrixDot->setBrush(QBrush(QColor(255, 255, 0, 255)));
				ADO.Scene->addItem(MatrixDot);
			}
        }
    }

    vector<QColor> rgb;
    rgb.push_back(QColor(255,255,0));
    rgb.push_back(QColor(0,255,255));
    rgb.push_back(QColor(100,255,0));

	for(int i = 0; i < accepted_orfs.size() && i < count; i++)
	{
		int orf_id = accepted_orfs[i];
		floParam flp;
        flp.M = translations[sorted_translations[orf_id][2][0]];
        flp.start_w_atg = false;
        flp.is_complete = true;
		auto ADO_orf = ADO;
		ADO_orf.mColor = QColor(0, 0, 0, 255);
        if(i < rgb.size())
        {
            ADO_orf.mColor = rgb[i];
        }
        ADO_orf.UpperCorner += QPointF(0, (i+1)*(i%2==1?1:-1)*2);
        Draw_ORF(flp, Reference_base, 0, ADO_orf);
	}
}


/*
* Libs compare drawer
*/

void Draw_ORFs_Compare_Libs(vector<AlignedRow>& matrix,
	vector<AlignedRow>& matrix2,
	floParam& MainPath,
	vector<TlessDNA>& Reference_base,
	size_t Reference_id,
	AbstractDrawOptions ADO,
    TAlignerOptions Program_Options)
{
    auto& rtless = Reference_base[Reference_id];
	map<int, vector<int> > sorted_alignments, sorted_alignments2;

	for(size_t i = 0; i < matrix.size(); i++)
		if(matrix[i].ref == Reference_id)
			sorted_alignments[matrix[i].rf_start].push_back(i);

	for(size_t i = 0; i < matrix2.size(); i++)
		if(matrix2[i].ref == Reference_id)
			sorted_alignments2[matrix2[i].rf_start].push_back(i);

	map<int, map<int, int> > Editing_Matrix, Editing_Matrix2;
	map<int, int> Editing_Matrix_CC, Editing_Matrix_CC2;
	for(size_t i = 0; i < rtless.T.size(); i++)
		if(sorted_alignments.count(i))
			for(auto& id : sorted_alignments[i]) {
				auto& alignment = matrix[id];
				for(size_t p = 0; p < alignment.r.size(); p++)
				{
					Editing_Matrix[alignment.r[p]][i+p]++;
					Editing_Matrix_CC[i+p]++;
				}

			}
	for(size_t i = 0; i < rtless.T.size(); i++)
		if(sorted_alignments2.count(i))
			for(auto& id : sorted_alignments2[i]) {
				auto& alignment = matrix2[id];
				for(size_t p = 0; p < alignment.r.size(); p++)
				{
					Editing_Matrix2[alignment.r[p]][i+p]++;
					Editing_Matrix_CC2[i+p]++;
				}
			}

	map<int, map<int, int> > MainXY;
	taORF CDS = Find_Longest_ORF(MainPath);
	int MainXS = MainPath.M.rf_start + CDS.tlstart;
	int MainXE = MainXS + CDS.tlend-CDS.tlstart;
	for(int i = 0; i < MainPath.M.mp.dT.size(); i++)
	{
		MainXY[MainPath.M.rf_start+i][MainPath.M.mp.dT[i+1] - rtless.dT[MainPath.M.rf_start+i+1]] = 1;
	}

	int counter_dots_diff_red = 0;
	int counter_dots_diff_blue = 0;
	int counter_dots_diff_red_ref = 0;
	int counter_dots_diff_blue_ref = 0;
	int counter_dots_diff_red_main = 0;
	int counter_dots_diff_blue_main = 0;
	int counter_dots_diff_red_alt = 0;
	int counter_dots_diff_blue_alt = 0;


	map<int, map<int, QColor> > Actual_dot_alpha;
	double norm_coeff = 1.0 * matrix.size() / matrix2.size();

	for(int i = 0; i < rtless.T.size(); i++)
		for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
			double dist = (1.0*Editing_Matrix[k][i]) / (norm_coeff*Editing_Matrix2[k][i]);
			//double dist = (1.0*Editing_Matrix[k][i]/Editing_Matrix_CC[i]) / (1.0*Editing_Matrix2[k][i]/Editing_Matrix_CC2[i]);
			if(Editing_Matrix[k][i] <= 4 || Editing_Matrix2[k][i] <= 4) dist = 1.0;
			if(dist < 3.0000 && dist > 0.3333)
				Actual_dot_alpha[k][i] = QColor(95, 105, 105);
			else if(dist >= 3.0000)
			{
				Actual_dot_alpha[k][i] = QColor(255, 0, 15, 100);
				if(i >= MainXS && i <= MainXE)
				{
					counter_dots_diff_red++;
					if(MainXY[i][k] == 1 && k != 0)
						counter_dots_diff_red_main++;
					else if(k == 0 && MainXY[i][k] != 1)
						counter_dots_diff_red_ref++;
					else if(k == 0 && MainXY[i][k] == 1) {}
					else counter_dots_diff_red_alt++;
				}
			}
			else if(dist <= 0.3333)
			{
				Actual_dot_alpha[k][i] = QColor(0, 15, 255, 100);
				if(i >= MainXS && i <= MainXE)
				{
					counter_dots_diff_blue++;
					if(MainXY[i][k] == 1 && k != 0)
						counter_dots_diff_blue_main++;
					else if(k == 0 && MainXY[i][k] != 1)
						counter_dots_diff_blue_ref++;
					else if(k == 0 && MainXY[i][k] == 1) {}
					else counter_dots_diff_blue_alt++;
				}
			}
		}

	cout << "\n--- Compare stats ---\n\tRed:\n\t\tAll\t" <<
		counter_dots_diff_red << "\n\t\tRef\t" << counter_dots_diff_red_ref << "\n\t\tMain\t" <<
		counter_dots_diff_red_main << "\n\t\tAlt\t" << counter_dots_diff_red_alt << "\n\tBlue:\n\t\tAll\t" <<
		counter_dots_diff_blue << "\n\t\tRef\t" << counter_dots_diff_blue_ref << "\n\t\tMain\t" <<
		counter_dots_diff_blue_main << "\n\t\tAlt\t" << counter_dots_diff_blue_alt << "\n\n";

    // Draw matrix cycle
    double pos_x = 0;
    for(int i = 0; i < rtless.T.size(); i++)
    {
        pos_x += x_axis_step;
        for(int k = -max_indel_l; k <= max_indel_l; k++)
        {
            int cell_rd_count = Editing_Matrix[k][i];
            if(cell_rd_count >= draw_dot_rdsupp_min)
            {
				QGraphicsEllipseItem* MatrixDot = new QGraphicsEllipseItem(
					pos_x + 1.0,
					ADO.UpperCorner.y() + x_axis_step*(-1)*k - 1.5*x_axis_step + 1.0, // - 2.0,
					6.0, 6.0);
				MatrixDot->setPen(QPen(Actual_dot_alpha[k][i]));
				MatrixDot->setBrush(QBrush(Actual_dot_alpha[k][i]));
				ADO.Scene->addItem(MatrixDot);
            }
        }
    }
}



vector<int> Recalc_Support(MappedPart& ORF,
    vector<MappedPart>& Alignments, int& supp_total, int& supp_edited)
{
    vector<int> selected_rd;
	supp_total = 0; supp_edited = 0;
    for(int i = 0; i < Alignments.size(); i++)
    {
        auto& read = Alignments[i];
        if(read.rf_start >= ORF.rf_start) {
            bool exact = true;
            int delta = read.rf_start - ORF.rf_start;
            for(int j = 1; j < ORF.mp.dT.size()-1 && j < read.mp.dT.size()-1; j++)
            {
                if(ORF.mp.dT[j+delta] != read.mp.dT[j]) {
                    exact = false; break;
                }
            }

            if(exact) {
                selected_rd.push_back(i);
				supp_total++;
				if(!read.as_ref) supp_edited++;
            }
        }
    }

    return selected_rd;
}

void Draw_ORF_Mapped_Reads(MappedPart& ORF,
    vector<MappedPart>& Alignments,
	vector<TlessDNA>& Reference_base,
	size_t Reference_id,
	AbstractDrawOptions ADO)
{
	int edited_supp, total_supp;
    auto AlignedIDs = Recalc_Support(ORF, Alignments, total_supp, edited_supp);

    static double curr_y = ADO.UpperCorner.y();
	int curr_y_begin = curr_y;
    for(int i = 0; i < AlignedIDs.size(); i++)
    {
        int id = AlignedIDs[i];
		if(Alignments[id].as_ref) continue;
        QGraphicsRectItem* rdl = new QGraphicsRectItem(x_axis_step*Alignments[id].rf_start, (++curr_y)++,
            x_axis_step*Alignments[id].mp.T.size(), 1);
        rdl->setPen(QPen(ADO.mColor));
        rdl->setBrush(QBrush(ADO.mColor));
        ADO.Scene->addItem(rdl);
    }
	/*
		Calculate site coverage
	*/
	map<int, int> coverage_vector;
	for(int i = 0; i < AlignedIDs.size(); i++)
	{
		for(int j = 0; j < Alignments[AlignedIDs[i]].mp.T.size(); j++)
		{
			coverage_vector[Alignments[AlignedIDs[i]].rf_start+j]++;
		}
	}
	int badly_covered_sites = 0; for(auto& x : coverage_vector) if(x.second <= 1) badly_covered_sites++;
	/*
		Draw legend
	*/
	QFont Font1("Arial", 30, QFont::Bold);
	auto frY = curr_y_begin + (curr_y - curr_y_begin) / 3;
	//QGraphicsRectItem* borderRect = new QGraphicsRectItem(
	//	x_axis_step * (Reference_base[0].T.size() - 15),
	//	frY, 70, 30);
	//ADO.Scene->addItem(borderRect);
	auto text = ADO.Scene->addText("" +
		QString(i2s(total_supp).c_str()) + " | " +
		QString(i2s(edited_supp).c_str()) + " | " +
		QString(i2s(badly_covered_sites).c_str()), Font1);
	text->setDefaultTextColor(ADO.mColor);
    text->setPos(x_axis_step * (Reference_base[0].T.size() - 20), frY + 2);
}

void Draw_Alternatively_Edited_Sites(
	MappedPart& MP,
	MappedPart& XX,
	vector<TlessDNA>& Reference_base,
	size_t Reference_id,
	AbstractDrawOptions ADO)
{
	floParam flp1;
	flp1.M = MP;
	flp1.start_w_atg = false;
	flp1.is_complete = true;
	auto OMP = Find_Longest_ORF(flp1);
	floParam flp2;
	flp2.M = XX;
	flp2.start_w_atg = false;
	flp2.is_complete = true;
	auto OXX = Find_Longest_ORF(flp2);

	int xx_orf_rf_start = OXX.tlstart + XX.rf_start;
	int mp_orf_rf_start = OMP.tlstart + MP.rf_start;
	int xx_dt_shift, mp_dt_shift, rf_alignment_start;

	if(mp_orf_rf_start >= xx_orf_rf_start)
	{
		rf_alignment_start = mp_orf_rf_start;
		mp_dt_shift = OMP.tlstart;
		xx_dt_shift = rf_alignment_start - XX.rf_start;
	}
	else
	{
		rf_alignment_start = xx_orf_rf_start;
		xx_dt_shift = OXX.tlstart;
		mp_dt_shift = rf_alignment_start - MP.rf_start;
	}

	int debug_count_circles = 0;
	for(int i = 0; i < OXX.tlend-xx_dt_shift  /*OXX.tlstart?*/ && i < OMP.tlend-mp_dt_shift /*OMP.tlstart*/; i++)
	{
		if(XX.mp.dT[xx_dt_shift+i+1] != MP.mp.dT[mp_dt_shift+i+1])
		{
			QGraphicsEllipseItem* aes = new QGraphicsEllipseItem(
				x_axis_step*(1+rf_alignment_start+i),
				ADO.UpperCorner.y(),
				x_axis_step,
	            10);
	        aes->setPen(QPen(ADO.mColor));
	        aes->setBrush(QBrush(ADO.mColor));
	        ADO.Scene->addItem(aes);

			debug_count_circles++;
		}
	}

	//cout << "debug_count_circles = " << debug_count_circles << "\n";
}



int Get_Total_ES_dist(MappedPart& MP, MappedPart& XX,
	vector<TlessDNA>& Reference_base, size_t Reference_id,
	vector<MappedPart>& Alignments, double& ess_ratio, int& ess_supp)
{


	floParam flp1;
	flp1.M = MP;
	flp1.start_w_atg = false;
	flp1.is_complete = true;
	auto OMP = Find_Longest_ORF(flp1);
	floParam flp2;
	flp2.M = XX;
	flp2.start_w_atg = false;
	flp2.is_complete = true;
	auto OXX = Find_Longest_ORF(flp2);

	int xx_orf_rf_start = OXX.tlstart + XX.rf_start;
	int mp_orf_rf_start = OMP.tlstart + MP.rf_start;
	int xx_dt_shift, mp_dt_shift, rf_alignment_start;

	if(mp_orf_rf_start >= xx_orf_rf_start)
	{
		rf_alignment_start = mp_orf_rf_start;
		mp_dt_shift = OMP.tlstart;
		xx_dt_shift = rf_alignment_start - XX.rf_start;
	}
	else
	{
		rf_alignment_start = xx_orf_rf_start;
		xx_dt_shift = OXX.tlstart;
		mp_dt_shift = rf_alignment_start - MP.rf_start;
	}

	int ESD = 0;
	vector<vector<int> > ess;
	for(int i = 0; i < OXX.tlend-xx_dt_shift  /*OXX.tlstart?*/ && i < OMP.tlend-mp_dt_shift /*OMP.tlstart*/; i++)
	{
		if(XX.mp.dT[xx_dt_shift+i+1] != MP.mp.dT[mp_dt_shift+i+1])
		{
			ESD++;
			vector<int> d2;
			d2.push_back(xx_dt_shift+i+1);
			d2.push_back(mp_dt_shift+i+1);
			ess.push_back(d2);
		}
	}

	vector<vector<int> > supp;

	for(int i = 0; i < ess.size(); i++)
	{
		vector<int> s2(2,0);
		for(int a = 0; a < Alignments.size(); a++)
		{
			auto& read = Alignments[a];

			if(read.rf_start >= XX.rf_start) {
	            bool exact = true;
	            int delta = read.rf_start - XX.rf_start;
	            for(int j = 1; j < XX.mp.dT.size()-1 && j < read.mp.dT.size()-1; j++)
	            {
	                if(XX.mp.dT[j+delta] != read.mp.dT[j]) {
	                    exact = false; break;
	                }
	            }

	            if(exact && read.rf_start + read.mp.dT.size() > XX.rf_start + ess[i][0] && read.rf_start < XX.rf_start + ess[i][0]) {
					s2[0]++;
	            }
	        }

			if(read.rf_start >= MP.rf_start) {
	            bool exact = true;
	            int delta = read.rf_start - MP.rf_start;
	            for(int j = 1; j < MP.mp.dT.size()-1 && j < read.mp.dT.size()-1; j++)
	            {
	                if(MP.mp.dT[j+delta] != read.mp.dT[j]) {
	                    exact = false; break;
	                }
	            }

	            if(exact && read.rf_start + read.mp.dT.size() > MP.rf_start + ess[i][1] && read.rf_start < MP.rf_start + ess[i][1]) {
					s2[1]++;
	            }
	        }
		}

		supp.push_back(s2);
	}

	double mp_min = 9999999;
	double xx_min = 9999999;
	for(int i = 0; i < supp.size(); i++)
	{
		if(xx_min > supp[i][0]) xx_min = supp[i][0];
		if(mp_min > supp[i][1]) mp_min = supp[i][1];
	}

	if(supp.size() == 0) ess_supp = -1;
	ess_ratio = xx_min / mp_min;
	ess_supp = xx_min;
	return ESD;
}


int Get_Total_ES_dist_mRNA(MappedPart& MP, MappedPart& XX,
	vector<TlessDNA>& Reference_base, size_t Reference_id,
	vector<MappedPart>& Alignments, double& ess_ratio, int& ess_supp)
{


	int xx_dt_shift, mp_dt_shift, rf_alignment_start;

	if(MP.rf_start >= XX.rf_start)
	{
		rf_alignment_start = MP.rf_start;
		mp_dt_shift = 0;
		xx_dt_shift = rf_alignment_start - XX.rf_start;
	}
	else
	{
		rf_alignment_start = XX.rf_start;
		xx_dt_shift = 0;
		mp_dt_shift = rf_alignment_start - MP.rf_start;
	}

	int ESD = 0;
	vector<vector<int> > ess;
	for(int i = 1; i < XX.mp.dT.size() - xx_dt_shift - 1 && i < MP.mp.dT.size() - mp_dt_shift - 1; i++)
	{
		if(XX.mp.dT[xx_dt_shift+i+1] != MP.mp.dT[mp_dt_shift+i+1])
		{
			ESD++;
			vector<int> d2;
			d2.push_back(xx_dt_shift+i+1);
			d2.push_back(mp_dt_shift+i+1);
			ess.push_back(d2);
		}
	}

	vector<vector<int> > supp;

	for(int i = 0; i < ess.size(); i++)
	{
		vector<int> s2(2,0);
		for(int a = 0; a < Alignments.size(); a++)
		{
			auto& read = Alignments[a];

			if(read.rf_start >= XX.rf_start) {
	            bool exact = true;
	            int delta = read.rf_start - XX.rf_start;
	            for(int j = 1; j < XX.mp.dT.size()-3 && j < read.mp.dT.size()-3; j++)
	            {
	                if(XX.mp.dT[j+delta] != read.mp.dT[j]) {
	                    exact = false; break;
	                }
	            }

	            if(exact && read.rf_start + read.mp.dT.size() > XX.rf_start + ess[i][0] && read.rf_start < XX.rf_start + ess[i][0]) {
					s2[0]++;
	            }
	        }

			if(read.rf_start >= MP.rf_start) {
	            bool exact = true;
	            int delta = read.rf_start - MP.rf_start;
	            for(int j = 1; j < MP.mp.dT.size()-3 && j < read.mp.dT.size()-3; j++)
	            {
	                if(MP.mp.dT[j+delta] != read.mp.dT[j]) {
	                    exact = false; break;
	                }
	            }

	            if(exact && read.rf_start + read.mp.dT.size() > MP.rf_start + ess[i][1] && read.rf_start < MP.rf_start + ess[i][1]) {
					s2[1]++;
	            }
	        }
		}

		supp.push_back(s2);
	}

	double mp_min = 9999999;
	double xx_min = 9999999;
	for(int i = 0; i < supp.size(); i++)
	{
		if(xx_min > supp[i][0]) xx_min = supp[i][0];
		if(mp_min > supp[i][1]) mp_min = supp[i][1];
	}

	if(supp.size() == 0) ess_supp = -1;
	ess_ratio = xx_min / mp_min;
	ess_supp = xx_min;
	return ESD;
}






int Count_Alt_Edited_Sites(MappedPart& XX, MappedPart& MP, TlessDNA& RF)
{
	int xx_dt_shift, mp_dt_shift, rf_alignment_start;

	if(MP.rf_start >= XX.rf_start)
	{
		if(XX.rf_start+XX.mp.dT.size()-1 < MP.rf_start) return 0;
		rf_alignment_start = MP.rf_start;
		mp_dt_shift = 0;
		xx_dt_shift = rf_alignment_start - XX.rf_start;
	}
	else
	{
		if(MP.rf_start+MP.mp.dT.size()-1 < XX.rf_start) return 0;
		rf_alignment_start = XX.rf_start;
		xx_dt_shift = 0;
		mp_dt_shift = rf_alignment_start - MP.rf_start;
	}

	int ESD = 0;
	for(int i = 0; i < XX.mp.dT.size() - xx_dt_shift - 2 && i < MP.mp.dT.size() - mp_dt_shift - 2; i++)
	{
		if(XX.mp.dT[xx_dt_shift+i+1] != MP.mp.dT[mp_dt_shift+i+1])
		{
			ESD++;
			if(XX.mp.dT[xx_dt_shift+i+1]-RF.dT[XX.rf_start+xx_dt_shift+i+1] == 0)
			{
				ESD--;
			}
		}
	}
	return ESD;
}



#endif
