#include <fstream>
#include <iomanip>
#include <iostream>

#include <QGridLayout>
#include <QHBoxLayout>
#include <QMessageBox>
#include <QVBoxLayout>

#include "outWindow.hpp"

/*------------------------------------------------------------------------------
 * Constructor of the class OutWindow.
 *
 * Inputs:
 *  - nFOutData: directory containing the output files.
------------------------------------------------------------------------------*/

OutWindow::OutWindow(std::string nFOutData) : QWidget(){
    m_white      = QColor(255, 255, 255);
    m_blueTum    = QColor(46, 165, 225);
    m_red        = QColor(192, 0, 0);
    m_redTum     = QColor(237, 125, 49);
    m_black      = QColor(56, 56, 56);
    m_blueTumDam = QColor(0, 32, 96);
    m_blueG1     = QColor(46, 165, 225);
    m_green      = QColor(156, 204, 88);
    m_orange     = QColor(246, 161, 27);
    m_violet     = QColor(109, 96, 214);
    m_redG0      = QColor(187, 81, 51);
    m_yellowNec  = QColor(127, 96, 0);
    m_brownMit   = QColor(132, 60, 12);
    m_greenApop  = QColor(56, 87, 35);

    m_selChartGroup = new QGroupBox("Select a chart", this);
    m_selMapGroup   = new QGroupBox("Select a map", this);
    m_mapGroup      = new QGroupBox(this);

    m_selChart = new QComboBox(m_selChartGroup);
    m_selMap   = new QComboBox(m_selMapGroup);

    m_chartView = new QChartView;
    m_sDash = new QLineSeries;
    m_sDash->append(0.0, 0.0);
    m_sDash->append(0.0, 100.0);
    m_sDash->setPen(QPen(Qt::DashLine));

    std::ifstream fEndTreatTime((nFOutData + "/endTreatTumDens.res").c_str());

    if(!fEndTreatTime.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening "
                                             "endTreatTumDens.res");
    }

    else{
        fEndTreatTime >> m_endTreatTime;
        fEndTreatTime.close();
    }



    std::ifstream fRec((nFOutData + "/rec.res").c_str());

    if(!fRec.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening rec.res");
    }

    else{
        fRec >> m_rec;
        if (m_rec){
            fRec >> m_recTime;
            fRec.close();
        }
    }

    m_endTreatDash = new QLineSeries;
    m_endTreatDash->append(m_endTreatTime, 0.0);
    m_endTreatDash->append(m_endTreatTime, 100.0);
    QPen pen;
    pen.setStyle(Qt::DashLine);
    pen.setColor(m_green);
    m_endTreatDash->setPen(pen);
    if (m_rec){
        m_recDash = new QLineSeries;
        m_recDash->append(m_recTime, 0.0);
        m_recDash->append(m_recTime, 100.0);
        pen.setColor(m_red);
        m_recDash->setPen(pen);
    }

    double a, b, c, d, e, f;

    std::ifstream fTumDens((nFOutData + "/tumDens.res").c_str());

    if(!fTumDens.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening "
                                             "tumDens.res");
    }

    else{
        m_cTumDens = new QChart;
        m_xTumDens = new QValueAxis;
        m_yTumDens = new QValueAxis;
        m_cTumDens->addAxis(m_xTumDens, Qt::AlignBottom);
        m_cTumDens->addAxis(m_yTumDens, Qt::AlignLeft);
        m_cTumDens->setMaximumWidth(1000);

        QLineSeries *sTumDens = new QLineSeries(m_cTumDens);

        fTumDens >> a >> b;
        while(!fTumDens.eof()){
            sTumDens->append(a, b);
            fTumDens >> a >> b;
        }
        fTumDens.close();
        QPen pen;
        pen.setWidth(5);
        pen.setColor(m_blueTum);
        sTumDens->setPen(pen);
        m_cTumDens->addSeries(sTumDens);
        m_cTumDens->legend()->hide();
        m_cTumDens->setTitle("Evolution of the tumor density");

        m_xTumDens->setTitleText("Time (h)");
        m_xTumDens->setLabelFormat("%i");
        sTumDens->attachAxis(m_xTumDens);
        m_xTumDens->applyNiceNumbers();

        m_yTumDens->setTitleText("Tumor density");
        m_yTumDens->setLabelFormat("%i");
        sTumDens->attachAxis(m_yTumDens);
        m_yTumDens->setMin(0.0);
        m_yTumDens->applyNiceNumbers();

        m_cTumDens->addSeries(m_sDash);
        m_sDash->attachAxis(m_xTumDens);
        m_sDash->attachAxis(m_yTumDens);

        m_cTumDens->addSeries(m_endTreatDash);
        m_endTreatDash->attachAxis(m_xTumDens);
        m_endTreatDash->attachAxis(m_yTumDens);

        if(m_rec){
            m_cTumDens->addSeries(m_recDash);
            m_recDash->attachAxis(m_xTumDens);
            m_recDash->attachAxis(m_yTumDens);
        }

        m_selChart->addItem("Tumor density");
    }

    std::ifstream fTumVol((nFOutData + "/tumVol.res").c_str());

    if(!fTumVol.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening "
                                             "tumVol.res");
    }

    else{
        m_cTumVol = new QChart;
        m_xTumVol = new QValueAxis;
        m_yTumVol = new QValueAxis;
        m_cTumVol->addAxis(m_xTumVol, Qt::AlignBottom);
        m_cTumVol->addAxis(m_yTumVol, Qt::AlignLeft);
        m_cTumVol->setMaximumWidth(1000);

        QLineSeries *sTumVol = new QLineSeries(m_cTumVol);

        fTumVol >> a >> b;
        while(!fTumVol.eof()){
            sTumVol->append(a, b);
            fTumVol >> a >> b;
        }
        fTumVol.close();

        m_cTumVol->addSeries(sTumVol);
        m_cTumVol->legend()->hide();
        m_cTumVol->setTitle("Evolution of the tumor volume");

        m_xTumVol->setTitleText("Time (h)");
        m_xTumVol->setLabelFormat("%i");
        sTumVol->attachAxis(m_xTumVol);
        m_xTumVol->applyNiceNumbers();

        m_yTumVol->setTitleText("Tumor volume (mm³)");
        m_yTumVol->setLabelFormat("%4.2f");
        sTumVol->attachAxis(m_yTumVol);
        m_yTumVol->setMin(0.0);
        m_yTumVol->applyNiceNumbers();

        m_selChart->addItem("Tumor volume");
    }

    std::ifstream fVascDens((nFOutData + "/vascDens.res").c_str());

    if(!fVascDens.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening "
                                             "vascDens.res");
    }

    else{
        m_cVascDens = new QChart;
        m_xVascDens = new QValueAxis;
        m_yVascDens = new QValueAxis;
        m_cVascDens->addAxis(m_xVascDens, Qt::AlignBottom);
        m_cVascDens->addAxis(m_yVascDens, Qt::AlignLeft);
        m_cVascDens->setMaximumWidth(1000);

        QLineSeries *sVascDens     = new QLineSeries(m_cVascDens);
        QLineSeries *sNormVascDens = new QLineSeries(m_cVascDens);
        QLineSeries *sTumVascDens  = new QLineSeries(m_cVascDens);


        fVascDens >> a >> b >> c >> d;
        while(!fVascDens.eof()){
            sVascDens->append(a, b);
            sNormVascDens->append(a, c);
            sTumVascDens->append(a, d);
            fVascDens >> a >> b >> c >> d;
        }
        fVascDens.close();

        m_cVascDens->addSeries(sVascDens);
        QPen pen(m_red);
        pen.setWidth(2);
        sNormVascDens->setPen(pen);
        m_cVascDens->addSeries(sNormVascDens);
        pen.setColor(m_redTum);
        sTumVascDens->setPen(pen);
        m_cVascDens->addSeries(sTumVascDens);

        m_cVascDens->setTitle("Evolution of the vascular density");
        sVascDens->setName("Total");
        sNormVascDens->setName("Normal vessels");
        sTumVascDens->setName("Tumor vessels");

        m_xVascDens->setTitleText("Time (h)");
        m_xVascDens->setLabelFormat("%i");
        sVascDens->attachAxis(m_xVascDens);
        sNormVascDens->attachAxis(m_xVascDens);
        sTumVascDens->attachAxis(m_xVascDens);
        m_xVascDens->applyNiceNumbers();

        m_yVascDens->setTitleText("Vascular density");
        m_yVascDens->setLabelFormat("%4.2f");
        sVascDens->attachAxis(m_yVascDens);
        sNormVascDens->attachAxis(m_yVascDens);
        sTumVascDens->attachAxis(m_yVascDens);
        m_yVascDens->setMin(0.0);
        m_yVascDens->applyNiceNumbers();

        m_selChart->addItem("Vascular density");
    }

    std::ifstream fKilledCells((nFOutData + "/killedCells.res").c_str());

    if(!fKilledCells.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening "
                                             "killedCells.res");
    }

    else{
        m_cKilledCells = new QChart;
        m_xKilledCells = new QValueAxis;
        m_yKilledCells = new QValueAxis;
        m_cKilledCells->addAxis(m_xKilledCells, Qt::AlignBottom);
        m_cKilledCells->addAxis(m_yKilledCells, Qt::AlignLeft);
        m_cKilledCells->setMaximumWidth(1000);

        QLineSeries *sKilledCells = new QLineSeries(m_cKilledCells);

        fKilledCells >> a >> b;
        while(!fKilledCells.eof()){
            sKilledCells->append(a, b);
            fKilledCells >> a >> b;
        }
        fKilledCells.close();

        m_cKilledCells->addSeries(sKilledCells);
        m_cKilledCells->legend()->hide();
        m_cKilledCells->setTitle("Evolution of the percentage of killed tumor "
                                 "cells");

        m_xKilledCells->setTitleText("Time (h)");
        m_yKilledCells->setLabelFormat("%i");
        sKilledCells->attachAxis(m_xKilledCells);
        m_xKilledCells->applyNiceNumbers();

        m_yKilledCells->setTitleText("Percentage of killed tumor cells");
        m_yKilledCells->setLabelFormat("%i");
        sKilledCells->attachAxis(m_yKilledCells);
        m_yKilledCells->setMin(0.0);
        m_yKilledCells->applyNiceNumbers();

        m_selChart->addItem("Percentage of killed tumor cells");
    }

    std::ifstream fDeadDens((nFOutData + "/deadCellsDens.res").c_str());

    if(!fDeadDens.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening "
                                             "deadCellsDens.res");
    }

    else{
        m_cDeadDens = new QChart;
        m_xDeadDens = new QValueAxis;
        m_yDeadDens = new QValueAxis;
        m_cDeadDens->addAxis(m_xDeadDens, Qt::AlignBottom);
        m_cDeadDens->addAxis(m_yDeadDens, Qt::AlignLeft);
        m_cDeadDens->setMaximumWidth(1000);

        QLineSeries *sDeadDens = new QLineSeries(m_cDeadDens);

        fDeadDens >> a >> b;
        while(!fDeadDens.eof()){
            sDeadDens->append(a, b);
            fDeadDens >> a >> b;
        }
        fDeadDens.close();

        m_cDeadDens->addSeries(sDeadDens);
        m_cDeadDens->legend()->hide();
        m_cDeadDens->setTitle("Evolution of the dead cell density");

        m_xDeadDens->setTitleText("Time (h)");
        m_xDeadDens->setLabelFormat("%i");
        sDeadDens->attachAxis(m_xDeadDens);
        m_xDeadDens->applyNiceNumbers();

        m_yDeadDens->setTitleText("Dead cell density");
        m_yDeadDens->setLabelFormat("%4.2f");
        sDeadDens->attachAxis(m_yDeadDens);
        m_yDeadDens->setMin(0.0);
        m_yDeadDens->applyNiceNumbers();

        m_selChart->addItem("Dead cell density");
    }

    std::ifstream fCycle((nFOutData + "/cycle.res").c_str());

    if(!fCycle.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening cycle.res");
    }

    else{
        m_cCycle = new QChart;
        m_xCycle = new QValueAxis;
        m_yCycle = new QValueAxis;
        m_cCycle->addAxis(m_xCycle, Qt::AlignBottom);
        m_cCycle->addAxis(m_yCycle, Qt::AlignLeft);
        m_cCycle->setMaximumWidth(1000);

        QLineSeries *sG1 = new QLineSeries();
        QLineSeries *sS  = new QLineSeries();
        QLineSeries *sG2 = new QLineSeries();
        QLineSeries *sM  = new QLineSeries();
        QLineSeries *sG0 = new QLineSeries();

        fCycle >> a >> b >> c >> d >> e >> f;
        while(!fCycle.eof()){
            sG1->append(a, b);
            sS->append(a, c);
            sG2->append(a, d);
            sM->append(a, e);
            sG0->append(a, f);
            fCycle >> a >> b >> c >> d >> e >> f;
        }
        fCycle.close();

        m_cCycle->addSeries(sG1);
        m_cCycle->addSeries(sS);
        m_cCycle->addSeries(sG2);
        m_cCycle->addSeries(sM);
        m_cCycle->addSeries(sG0);

        m_cCycle->setTitle("Evolution of the distribution of tumor cells");
        sG1->setName("G1");
        sS->setName("S");
        sG2->setName("G2");
        sM->setName("M");
        sG0->setName("G0");

        m_xCycle->setTitleText("Time(h)");
        m_xCycle->setLabelFormat("%i");
        sG1->attachAxis(m_xCycle);
        sS->attachAxis(m_xCycle);
        sG2->attachAxis(m_xCycle);
        sM->attachAxis(m_xCycle);
        sG0->attachAxis(m_xCycle);
        m_xCycle->applyNiceNumbers();

        m_yCycle->setTitleText("Percentage of tumor cells in each phase");
        m_yCycle->setLabelFormat("%i");
        sG1->attachAxis(m_yCycle);
        sS->attachAxis(m_yCycle);
        sG2->attachAxis(m_yCycle);
        sM->attachAxis(m_yCycle);
        sG0->attachAxis(m_yCycle);
        m_yCycle->setRange(0.0, 100.0);
        m_yCycle->applyNiceNumbers();

        m_selChart->addItem("Cell cycle distribution");
    }

    std::ifstream fHypDens((nFOutData + "/hypDens.res").c_str());

    if(!fHypDens.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening "
                                             "hypDens.res");
    }

    else{
        m_cHypDens = new QChart;
        m_xHypDens = new QValueAxis;
        m_yHypDens = new QValueAxis;
        m_cHypDens->addAxis(m_xHypDens, Qt::AlignBottom);
        m_cHypDens->addAxis(m_yHypDens, Qt::AlignLeft);
        m_cHypDens->setMaximumWidth(1000);

        QLineSeries *sHypDens = new QLineSeries(m_cHypDens);

        fHypDens >> a >> b;
        while(!fHypDens.eof()){
            sHypDens->append(a, b);
            fHypDens >> a >> b;
        }
        fHypDens.close();

        m_cHypDens->addSeries(sHypDens);
        m_cHypDens->legend()->hide();
        m_cHypDens->setTitle("Evolution of the hypoxic density");

        m_xHypDens->setTitleText("Time (h)");
        m_xHypDens->setLabelFormat("%i");
        sHypDens->attachAxis(m_xHypDens);
        m_xHypDens->applyNiceNumbers();

        m_yHypDens->setTitleText("Hypoxic density");
        m_yHypDens->setLabelFormat("%i");
        sHypDens->attachAxis(m_yHypDens);
        m_yHypDens->setMin(0.0);
        m_yHypDens->applyNiceNumbers();

        m_selChart->addItem("Hypoxic density");
    }

    std::ifstream fPO2Stat((nFOutData + "/pO2Stat.res").c_str());

    if(!fPO2Stat.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening "
                                             "pO2Stat.res");
    }

    else{
        m_cPO2Stat = new QChart;
        m_xPO2Stat = new QValueAxis;
        m_yPO2Stat = new QValueAxis;
        m_cPO2Stat->addAxis(m_xPO2Stat, Qt::AlignBottom);
        m_cPO2Stat->addAxis(m_yPO2Stat, Qt::AlignLeft);
        m_cPO2Stat->setMaximumWidth(1000);

        QLineSeries *sPO2Med  = new QLineSeries(m_cPO2Stat);
        QLineSeries *sPO2Mean = new QLineSeries(m_cPO2Stat);

        fPO2Stat >> a >> b >> c;
        while(!fPO2Stat.eof()){
            sPO2Med->append(a, b);
            sPO2Mean->append(a, c);
            fPO2Stat >> a >> b >> c;
        }
        fPO2Stat.close();

        m_cPO2Stat->addSeries(sPO2Med);
        m_cPO2Stat->addSeries(sPO2Mean);
        m_cPO2Stat->setTitle("Evolution of the pO2 statistics");
        sPO2Med->setName("Median");
        sPO2Mean->setName("Mean");

        m_xPO2Stat->setTitleText("Time (h)");
        m_xPO2Stat->setLabelFormat("%i");
        sPO2Mean->attachAxis(m_xPO2Stat);
        sPO2Med->attachAxis(m_xPO2Stat);
        m_xPO2Stat->applyNiceNumbers();

        m_yPO2Stat->setTitleText("Value");
        m_yPO2Stat->setLabelFormat("%4.2f");
        sPO2Mean->attachAxis(m_yPO2Stat);
        sPO2Med->attachAxis(m_yPO2Stat);
        m_yPO2Stat->setMin(0.0);
        m_yPO2Stat->applyNiceNumbers();

        m_selChart->addItem("pO2 statistics");
    }

    std::ifstream fVegfStat((nFOutData + "/vegfStat.res").c_str());

    if(!fVegfStat.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening "
                                             "vegfStat.res");
    }

    else{
        m_cVegfStat = new QChart;
        m_xVegfStat = new QValueAxis;
        m_yVegfStat = new QValueAxis;
        m_cVegfStat->addAxis(m_xVegfStat, Qt::AlignBottom);
        m_cVegfStat->addAxis(m_yVegfStat, Qt::AlignLeft);
        m_cVegfStat->setMaximumWidth(1000);

        QLineSeries *sVegfMed  = new QLineSeries(m_cVegfStat);
        QLineSeries *sVegfMean = new QLineSeries(m_cVegfStat);

        fVegfStat >> a >> b >> c;
        while(!fVegfStat.eof()){
            sVegfMed->append(a, b);
            sVegfMean->append(a, c);
            fVegfStat >> a >> b >> c;
        }
        fVegfStat.close();

        m_cVegfStat->addSeries(sVegfMed);
        m_cVegfStat->addSeries(sVegfMean);
        m_cVegfStat->setTitle("Evolution of the VEGF statistics");
        sVegfMed->setName("Median");
        sVegfMean->setName("Mean");

        m_xVegfStat->setTitleText("Time (h)");
        m_xVegfStat->setLabelFormat("%i");
        sVegfMean->attachAxis(m_xVegfStat);
        sVegfMed->attachAxis(m_xVegfStat);
        m_xVegfStat->applyNiceNumbers();

        m_yVegfStat->setTitleText("Value");
        m_yVegfStat->setLabelFormat("%4.2f");
        sVegfMean->attachAxis(m_yVegfStat);
        sVegfMed->attachAxis(m_yVegfStat);
        m_yVegfStat->setMin(0.0);
        m_yVegfStat->applyNiceNumbers();

        m_selChart->addItem("VEGF statistics");
    }

    m_chartView->setChart(m_cTumDens);

    m_map = new QLabel(m_mapGroup);

    std::ifstream fTissueDim((nFOutData + "/tissueDim.dat").c_str());

    if(!fTissueDim.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening "
                                             "tissueDim.dat");
    }

    else{
        fTissueDim >> m_mapNrow >> m_mapNcol;
        fTissueDim.close();
    }

    int nrowNcol(m_mapNrow * m_mapNcol);

    m_mapSclFac = 16;

    std::ifstream fSimParam((nFOutData + "/simParam.dat").c_str());

    if(!fSimParam.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening "
                                             "simParam.dat");
    }

    else{
        fSimParam >> m_simTime >> m_simTime >> m_simTimeStep;
        fSimParam.close();
    }


    std::ifstream fState((nFOutData + "/state.res").c_str());

    if(!fState.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening state.res");
    }

    else{
        int i(0), tempi;

        fState >> tempi;
        while(!fState.eof()){
            m_state.push_back(std::vector<int>());
            for(int j(0); j < nrowNcol; j++){
                m_state[i].push_back(tempi);
                fState >> tempi;
            }
            qApp->processEvents();
            i++;
        }
        fState.close();

        m_pixState = new QPixmap;
        m_legendSt = new QGroupBox(m_mapGroup);
        m_imState = new QImage(m_mapSclFac * m_mapNcol, m_mapSclFac * m_mapNrow,
                               QImage::Format_Indexed8);
        m_imState->setColor(0, m_white.rgb());
        m_imState->setColor(1, m_blueTum.rgb());
        m_imState->setColor(2, m_blueTumDam.rgb());
        m_imState->setColor(3, m_red.rgb());
        m_imState->setColor(4, m_redTum.rgb());
        m_imState->setColor(5, m_yellowNec.rgb());
        m_imState->setColor(6, m_brownMit.rgb());
        m_imState->setColor(7, m_greenApop.rgb());

        m_fib     = new QLabel("Healthy cells", m_legendSt);
        m_tum     = new QLabel("Tumor cells", m_legendSt);
        m_tumDam  = new QLabel("Damaged tumor cells", m_legendSt);
        m_normVes = new QLabel("Pre-existing vessels", m_legendSt);
        m_tumVes  = new QLabel("Neo-created vessels", m_legendSt);
        m_hypNec  = new QLabel("Hypoxic necrotic cells", m_legendSt);
        m_mitCat  = new QLabel("Mitotic catastrophe cells", m_legendSt);
        m_apop    = new QLabel("Apoptotic cells", m_legendSt);

        QImage fibIm     = QImage(m_mapSclFac, m_mapSclFac,
                                  QImage::Format_RGB32);
        QImage tumIm     = QImage(m_mapSclFac, m_mapSclFac,
                                  QImage::Format_RGB32);
        QImage tumDamIm  = QImage(m_mapSclFac, m_mapSclFac,
                                  QImage::Format_RGB32);
        QImage normVesIm = QImage(m_mapSclFac, m_mapSclFac,
                                  QImage::Format_RGB32);
        QImage tumVesIm  = QImage(m_mapSclFac, m_mapSclFac,
                                  QImage::Format_RGB32);
        QImage hypNecIm  = QImage(m_mapSclFac, m_mapSclFac,
                                  QImage::Format_RGB32);
        QImage mitCatIm  = QImage(m_mapSclFac, m_mapSclFac,
                                  QImage::Format_RGB32);
        QImage apopIm    = QImage(m_mapSclFac, m_mapSclFac,
                                  QImage::Format_RGB32);

        fibIm.fill(m_white.rgb());
        tumIm.fill(m_blueTum.rgb());
        tumDamIm.fill(m_blueTumDam.rgb());
        normVesIm.fill(m_red.rgb());
        tumVesIm.fill(m_redTum.rgb());
        hypNecIm.fill(m_yellowNec.rgb());
        mitCatIm.fill(m_brownMit.rgb());
        apopIm.fill(m_greenApop.rgb());

        m_fibSq     = new QLabel(m_legendSt);
        m_tumSq     = new QLabel(m_legendSt);
        m_tumDamSq  = new QLabel(m_legendSt);
        m_normVesSq = new QLabel(m_legendSt);
        m_tumVesSq  = new QLabel(m_legendSt);
        m_hypNecSq  = new QLabel(m_legendSt);
        m_mitCatSq  = new QLabel(m_legendSt);
        m_apopSq    = new QLabel(m_legendSt);

        m_fibSq->setPixmap(QPixmap::fromImage(fibIm));
        m_tumSq->setPixmap(QPixmap::fromImage(tumIm));
        m_tumDamSq->setPixmap(QPixmap::fromImage(tumDamIm));
        m_normVesSq->setPixmap(QPixmap::fromImage(normVesIm));
        m_tumVesSq->setPixmap(QPixmap::fromImage(tumVesIm));
        m_hypNecSq->setPixmap(QPixmap::fromImage(hypNecIm));
        m_mitCatSq->setPixmap(QPixmap::fromImage(mitCatIm));
        m_apopSq->setPixmap(QPixmap::fromImage(apopIm));

        QGridLayout *legStLayout = new QGridLayout();

        legStLayout->addWidget(m_fibSq, 0, 0);
        legStLayout->addWidget(m_tumSq, 1, 0);
        legStLayout->addWidget(m_tumDamSq, 2, 0);
        legStLayout->addWidget(m_normVesSq, 3, 0);
        legStLayout->addWidget(m_tumVesSq, 4, 0);
        legStLayout->addWidget(m_hypNecSq, 5, 0);
        legStLayout->addWidget(m_mitCatSq, 6, 0);
        legStLayout->addWidget(m_apopSq, 7, 0);
        legStLayout->addWidget(m_fib, 0, 1);
        legStLayout->addWidget(m_tum, 1, 1);
        legStLayout->addWidget(m_tumDam, 2, 1);
        legStLayout->addWidget(m_normVes, 3, 1);
        legStLayout->addWidget(m_tumVes, 4, 1);
        legStLayout->addWidget(m_hypNec, 5, 1);
        legStLayout->addWidget(m_mitCat, 6, 1);
        legStLayout->addWidget(m_apop, 7, 1);

        m_legendSt->setMinimumWidth(100);
        m_legendSt->setLayout(legStLayout);

        m_selMap->addItem("State of cells");
    }

    std::ifstream fTimer((nFOutData + "/timer.res").c_str());

    if(!fTimer.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening timer.res");
    }

    else{
        int i(0), tempi;

        fTimer >> tempi;
        while(!fTimer.eof()){
            m_timer.push_back(std::vector<int>());
            for(int j(0); j < nrowNcol; j++){
                m_timer[i].push_back(tempi);
                fTimer >> tempi;
            }
            qApp->processEvents();
            i++;
        }
        fTimer.close();

        m_pixTimer = new QPixmap;
        m_legendCyc = new QGroupBox(m_mapGroup);
        m_imTimer = new QImage(m_mapSclFac * m_mapNcol, m_mapSclFac * m_mapNrow,
                               QImage::Format_Indexed8);
        m_imTimer->setColor(0, m_white.rgb());
        m_imTimer->setColor(1, m_blueG1.rgb());
        m_imTimer->setColor(2, m_green.rgb());
        m_imTimer->setColor(3, m_orange.rgb());
        m_imTimer->setColor(4, m_violet.rgb());
        m_imTimer->setColor(5, m_redG0.rgb());

        m_G1 = new QLabel(QString::fromStdString("G1" + std::string(25, ' ')),
                          m_legendCyc);
        m_S  = new QLabel("S", m_legendCyc);
        m_G2 = new QLabel("G2", m_legendCyc);
        m_M  = new QLabel("M", m_legendCyc);
        m_G0 = new QLabel("G0", m_legendCyc);

        QImage G1Im = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage SIm  = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage G2Im = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage MIm  = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage G0Im = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);

        G1Im.fill(m_blueG1.rgb());
        SIm.fill(m_green.rgb());
        G2Im.fill(m_orange.rgb());
        MIm.fill(m_violet.rgb());
        G0Im.fill(m_redG0.rgb());

        m_G1Sq = new QLabel(m_legendCyc);
        m_SSq  = new QLabel(m_legendCyc);
        m_G2Sq = new QLabel(m_legendCyc);
        m_MSq  = new QLabel(m_legendCyc);
        m_G0Sq = new QLabel(m_legendCyc);

        m_G1Sq->setPixmap(QPixmap::fromImage(G1Im));
        m_SSq->setPixmap(QPixmap::fromImage(SIm));
        m_G2Sq->setPixmap(QPixmap::fromImage(G2Im));
        m_MSq->setPixmap(QPixmap::fromImage(MIm));
        m_G0Sq->setPixmap(QPixmap::fromImage(G0Im));

        QGridLayout *legCycLayout = new QGridLayout();
        legCycLayout->addWidget(m_G1Sq, 0, 0);
        legCycLayout->addWidget(m_SSq, 1, 0);
        legCycLayout->addWidget(m_G2Sq, 2, 0);
        legCycLayout->addWidget(m_MSq, 3, 0);
        legCycLayout->addWidget(m_G0Sq, 4, 0);
        legCycLayout->addWidget(m_G1, 0, 1);
        legCycLayout->addWidget(m_S, 1, 1);
        legCycLayout->addWidget(m_G2, 2, 1);
        legCycLayout->addWidget(m_M, 3, 1);
        legCycLayout->addWidget(m_G0, 4, 1);

        m_legendCyc->setMinimumWidth(100);
        m_legendCyc->setLayout(legCycLayout);

        m_selMap->addItem("Cell cycle");
    }

    std::ifstream fPO2((nFOutData + "/po2.res").c_str());

    if(!fPO2.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening po2.res");
    }

    else{
        int i(0);
        double tempd;

        m_maxpO2 = 0.0;

        fPO2 >> tempd;
        while(!fPO2.eof()){
            m_pO2.push_back(std::vector<double>());
            for(int j(0); j < nrowNcol; j++){
                m_pO2[i].push_back(tempd);
                if(tempd > m_maxpO2){
                    m_maxpO2 = tempd;
                }
                fPO2 >> tempd;
            }
            qApp->processEvents();
            i++;
        }
        fPO2.close();

        m_pixPO2 = new QPixmap;
        m_legendPO2 = new QGroupBox(m_mapGroup);
        m_imPO2 = new QImage(m_mapSclFac * m_mapNcol, m_mapSclFac * m_mapNrow,
                             QImage::Format_RGB32);
        m_pO2Bar = new QLabel(m_legendPO2);
        m_pO2ValMax = new QLabel(QString::number(m_maxpO2) + " mmHg");
        m_pO2ValMin = new QLabel("0 mmHg");

        const int nrowOxy(50), ncolOxy(720);
        int compH;
        QImage oxyIm = QImage(ncolOxy, nrowOxy, QImage::Format_RGB32);

        for(int i(0); i < nrowOxy; i++){
            for(int j(0); j < ncolOxy; j++){
                compH = (ncolOxy - j) / 3 % 241;
                oxyIm.setPixelColor(j, i, QColor::fromHsv(compH, 200, 255));
            }
        }
        m_pO2Bar->setPixmap(QPixmap::fromImage(oxyIm));

        QGridLayout *legOxyLayout = new QGridLayout();
        legOxyLayout->addWidget(m_pO2ValMin, 0, 0);
        legOxyLayout->addWidget(m_pO2Bar, 0, 1, 1, 2);
        legOxyLayout->addWidget(m_pO2ValMax, 0, 3);

        m_legendPO2->setMinimumWidth(100);
        m_legendPO2->setLayout(legOxyLayout);

        m_selMap->addItem("pO2");
    }

    std::ifstream fVegf((nFOutData + "/vegf.res").c_str());

    if(!fVegf.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening vegf.res");
    }

    else{
        int i(0);
        double tempd;

        m_maxvegf = 0.0;

        fVegf >> tempd;
        while(!fVegf.eof()){
            m_vegf.push_back(std::vector<double>());
            for(int j(0); j < nrowNcol; j++){
                m_vegf[i].push_back(tempd);
                if(tempd > m_maxvegf){
                    m_maxvegf = tempd;
                }
                fVegf >> tempd;
            }
            qApp->processEvents();
            i++;
        }
        fVegf.close();

        if(m_maxvegf == 0.0){
            m_maxvegf = 20.0;
        }

        m_pixVegf = new QPixmap;
        m_legendVegf = new QGroupBox(m_mapGroup);
        m_imVegf = new QImage(m_mapSclFac * m_mapNcol, m_mapSclFac * m_mapNrow,
                              QImage::Format_RGB32);
        m_vegfBar = new QLabel(m_legendVegf);
        m_vegfValMax = new QLabel(QString::number(m_maxvegf) + "mol/μm³");
        m_vegfValMin = new QLabel("0 mol/μm³");

        const int nrowVegf(50), ncolVegf(720);
        int compH;
        QImage vegfIm = QImage(ncolVegf, nrowVegf, QImage::Format_RGB32);

        for(int i(0); i < nrowVegf; i++){
            for(int j(0); j < ncolVegf; j++){
                compH = (ncolVegf - j) / 3 % 241;
                vegfIm.setPixelColor(j, i, QColor::fromHsv(compH, 200, 255));
            }
        }
        m_vegfBar->setPixmap(QPixmap::fromImage(vegfIm));

        QGridLayout *legVegfLayout = new QGridLayout();
        legVegfLayout->addWidget(m_vegfValMin, 0, 0);
        legVegfLayout->addWidget(m_vegfBar, 0, 1, 1, 2);
        legVegfLayout->addWidget(m_vegfValMax, 0, 3);

        m_legendVegf->setMinimumWidth(100);
        m_legendVegf->setLayout(legVegfLayout);

        m_selMap->addItem("VEGF");
    }

    m_time   = new QLabel("Time (h)", m_mapGroup);
    m_timeS  = new QSpinBox(m_mapGroup);
    m_play   = new QPushButton(m_mapGroup);
    m_slider = new QSlider(Qt::Horizontal, this);

    m_timeS->setMinimum(0);
    m_timeS->setMaximum((m_state.size() - 1) * m_simTimeStep);
    m_slider->setMinimum(0);
    m_slider->setMaximum((m_state.size() - 1) * m_simTimeStep);
    m_play->setIcon(QIcon("../Figures/play.png"));
    m_play->setFixedSize(40, 40);
    m_play->setIconSize(QSize(30, 30));
    m_stPlay = false;

    m_sel          = new QPushButton("Select a new mode", this);
    m_change       = new QPushButton("Change input and parameters", this);
    m_saveOutFiles = new QPushButton("Save output files", this);
    m_saveChart    = new QPushButton("Save chart", this);
    m_saveMap      = new QPushButton("Save map of current iteration", this);
    m_saveAllMaps  = new QPushButton("Save maps of every iteration", this);
    m_newSim       = new QPushButton("New simulation", this);

    QVBoxLayout *chartLayout = new QVBoxLayout;
    chartLayout->addWidget(m_selChart);
    m_selChartGroup->setLayout(chartLayout);

    QVBoxLayout *selMapLayout = new QVBoxLayout;
    selMapLayout->addWidget(m_selMap);
    m_selMapGroup->setLayout(selMapLayout);

    QHBoxLayout *timeLayout = new QHBoxLayout;
    timeLayout->addWidget(m_time);
    timeLayout->addWidget(m_timeS);

    QGridLayout *mapLayout = new QGridLayout;
    mapLayout->addWidget(m_map, 0, 0, 1, 2);
    mapLayout->addWidget(m_legendSt, 1, 0, 1, 2);
    mapLayout->addWidget(m_legendCyc, 1, 0, 1, 2);
    mapLayout->addWidget(m_legendPO2, 1, 0, 1, 2);
    mapLayout->addWidget(m_legendVegf, 1, 0, 1, 2);
    mapLayout->addWidget(m_slider, 2, 0);
    mapLayout->addLayout(timeLayout, 3, 0, 1, 2);
    mapLayout->addWidget(m_play, 2, 1, 1, 1, Qt::AlignCenter);
    m_mapGroup->setLayout(mapLayout);

    QHBoxLayout *hLayout = new QHBoxLayout;
    hLayout->addWidget(m_sel);
    hLayout->addWidget(m_change);
    hLayout->addWidget(m_saveOutFiles);
    hLayout->addWidget(m_saveChart);
    hLayout->addWidget(m_saveMap);
    hLayout->addWidget(m_saveAllMaps);
    hLayout->addWidget(m_newSim);

    QGridLayout *gridLayout = new QGridLayout;
    gridLayout->addWidget(m_selChartGroup, 0, 0);
    gridLayout->addWidget(m_chartView, 1, 0);
    gridLayout->addWidget(m_selMapGroup, 0, 1);
    gridLayout->addWidget(m_mapGroup, 1, 1);

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addLayout(gridLayout);
    layout->addLayout(hLayout);
    setLayout(layout);

    QObject::connect(m_selMap, SIGNAL(currentIndexChanged(int)), this,
                     SLOT(changeNumMap(int)));
    QObject::connect(m_selChart, SIGNAL(currentIndexChanged(int)), this,
                     SLOT(changeChart(int)));
    QObject::connect(m_slider, SIGNAL(valueChanged(int)), m_timeS,
                     SLOT(setValue(int)));
    QObject::connect(m_timeS, SIGNAL(valueChanged(int)), m_slider,
                     SLOT(setValue(int)));
    QObject::connect(m_timeS, SIGNAL(valueChanged(int)), this,
                     SLOT(changeIter(int)));
    QObject::connect(this, SIGNAL(updateSlider(int)), m_slider,
                     SLOT(setValue(int)));
    QObject::connect(m_play, SIGNAL(clicked()), this, SLOT(play()));
    QObject::connect(m_sel, SIGNAL(clicked()), this, SLOT(sel()));
    QObject::connect(m_change, SIGNAL(clicked()), this, SLOT(change()));
    QObject::connect(m_newSim, SIGNAL(clicked()),this, SLOT(newSim()));
    QObject::connect(m_saveOutFiles, SIGNAL(clicked()), this,
                     SLOT(saveOutFiles()));
    QObject::connect(m_saveChart, SIGNAL(clicked()), this, SLOT(saveChart()));
    QObject::connect(m_saveMap, SIGNAL(clicked()), this, SLOT(saveMap()));
    QObject::connect(m_saveAllMaps, SIGNAL(clicked()), this,
                     SLOT(saveAllMaps()));

    drawMap(m_selMap->currentIndex(), m_slider->value());

    setWindowTitle("Radiotherapy Simulator");
    setWindowIcon(QIcon("../Figures/logo.png"));
    showMaximized();
}


/*------------------------------------------------------------------------------
 * This slot goes back to the InWindow to change the model and simulation
 * parameters.
------------------------------------------------------------------------------*/

void OutWindow::change(){
    new InWindow;
    close();
}


/*------------------------------------------------------------------------------
 * This slot changes the current chart.
 * Inputs:
 *  - numChart: chart to be visualised (0, tumour density; 1, tumour volume;
 *  2, vascular density; 3, killed cells; 4, dead cell density; 5, cell cycle
 *  distribution; 6, hypoxic density; 7, pO2 statistics; 8, VEGF statistics).
------------------------------------------------------------------------------*/

void OutWindow::changeChart(const int numChart){
    switch(numChart){
    case 0:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_cTumDens->addSeries(m_sDash);
        m_cTumDens->addSeries(m_endTreatDash);
        m_sDash->attachAxis(m_xTumDens);
        m_sDash->attachAxis(m_yTumDens);
        m_endTreatDash->attachAxis(m_xTumDens);
        m_endTreatDash->attachAxis(m_yTumDens);
        if(m_rec){
            m_chartView->chart()->removeSeries(m_recDash);
            m_cTumDens->addSeries(m_recDash);
            m_recDash->attachAxis(m_xTumDens);
            m_recDash->attachAxis(m_yTumDens);
        }
        m_chartView->setChart(m_cTumDens);
        break;
    }

    case 1:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_cTumVol->addSeries(m_sDash);
        m_cTumVol->addSeries(m_endTreatDash);
        m_sDash->attachAxis(m_xTumVol);
        m_sDash->attachAxis(m_yTumVol);
        m_endTreatDash->attachAxis(m_xTumVol);
        m_endTreatDash->attachAxis(m_yTumVol);
        if(m_rec){
            m_chartView->chart()->removeSeries(m_recDash);
            m_cTumVol->addSeries(m_recDash);
            m_recDash->attachAxis(m_xTumVol);
            m_recDash->attachAxis(m_yTumVol);
        }
        m_chartView->setChart(m_cTumVol);
        break;
    }

    case 2:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_cVascDens->addSeries(m_sDash);
        m_cVascDens->addSeries(m_endTreatDash);
        m_sDash->attachAxis(m_xVascDens);
        m_sDash->attachAxis(m_yVascDens);
        m_endTreatDash->attachAxis(m_xVascDens);
        m_endTreatDash->attachAxis(m_yVascDens);
        if(m_rec){
            m_chartView->chart()->removeSeries(m_recDash);
            m_cVascDens->addSeries(m_recDash);
            m_recDash->attachAxis(m_xVascDens);
            m_recDash->attachAxis(m_yVascDens);
        }
        m_chartView->setChart(m_cVascDens);
        m_chartView->chart()->legend()->markers(m_sDash).front()->
                setVisible(false);
        m_chartView->chart()->legend()->markers(m_endTreatDash).front()->
                setVisible(false);
        if(m_rec){
            m_chartView->chart()->legend()->markers(m_recDash).front()->
                    setVisible(false);
        }
        break;
    }

    case 3:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_cKilledCells->addSeries(m_sDash);
        m_cKilledCells->addSeries(m_endTreatDash);
        m_sDash->attachAxis(m_xKilledCells);
        m_sDash->attachAxis(m_yKilledCells);
        m_endTreatDash->attachAxis(m_xKilledCells);
        m_endTreatDash->attachAxis(m_yKilledCells);
        if(m_rec){
            m_chartView->chart()->removeSeries(m_recDash);
            m_cKilledCells->addSeries(m_recDash);
            m_recDash->attachAxis(m_xKilledCells);
            m_recDash->attachAxis(m_yKilledCells);
        }
        m_chartView->setChart(m_cKilledCells);
        break;
    }

    case 4:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_cDeadDens->addSeries(m_sDash);
        m_cDeadDens->addSeries(m_endTreatDash);
        m_sDash->attachAxis(m_xDeadDens);
        m_sDash->attachAxis(m_yDeadDens);
        m_endTreatDash->attachAxis(m_xDeadDens);
        m_endTreatDash->attachAxis(m_yDeadDens);
        if(m_rec){
            m_chartView->chart()->removeSeries(m_recDash);
            m_cDeadDens->addSeries(m_recDash);
            m_recDash->attachAxis(m_xDeadDens);
            m_recDash->attachAxis(m_yDeadDens);
        }
        m_chartView->setChart(m_cDeadDens);
        break;
    }

    case 5:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_cCycle->addSeries(m_sDash);
        m_cCycle->addSeries(m_endTreatDash);
        m_sDash->attachAxis(m_xCycle);
        m_sDash->attachAxis(m_yCycle);
        m_endTreatDash->attachAxis(m_xCycle);
        m_endTreatDash->attachAxis(m_yCycle);
        if(m_rec){
            m_chartView->chart()->removeSeries(m_recDash);
            m_cCycle->addSeries(m_recDash);
            m_recDash->attachAxis(m_xCycle);
            m_recDash->attachAxis(m_yCycle);
        }
        m_chartView->setChart(m_cCycle);
        m_chartView->chart()->legend()->markers(m_sDash).front()->
                setVisible(false);
        m_chartView->chart()->legend()->markers(m_endTreatDash).front()->
                setVisible(false);
        if(m_rec){
            m_chartView->chart()->legend()->markers(m_recDash).front()->
                    setVisible(false);
        }
        break;
    }

    case 6:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_cHypDens->addSeries(m_sDash);
        m_cHypDens->addSeries(m_endTreatDash);
        m_sDash->attachAxis(m_xHypDens);
        m_sDash->attachAxis(m_yHypDens);
        m_endTreatDash->attachAxis(m_xHypDens);
        m_endTreatDash->attachAxis(m_yHypDens);
        if(m_rec){
            m_chartView->chart()->removeSeries(m_recDash);
            m_cHypDens->addSeries(m_recDash);
            m_recDash->attachAxis(m_xHypDens);
            m_recDash->attachAxis(m_yHypDens);
        }
        m_chartView->setChart(m_cHypDens);
        break;
    }

    case 7:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_cPO2Stat->addSeries(m_sDash);
        m_cPO2Stat->addSeries(m_endTreatDash);
        m_sDash->attachAxis(m_xPO2Stat);
        m_sDash->attachAxis(m_yPO2Stat);
        m_endTreatDash->attachAxis(m_xPO2Stat);
        m_endTreatDash->attachAxis(m_yPO2Stat);
        if(m_rec){
            m_chartView->chart()->removeSeries(m_recDash);
            m_cPO2Stat->addSeries(m_recDash);
            m_recDash->attachAxis(m_xPO2Stat);
            m_recDash->attachAxis(m_yPO2Stat);
        }
        m_chartView->setChart(m_cPO2Stat);
        m_chartView->chart()->legend()->markers(m_sDash).front()->
                setVisible(false);
        m_chartView->chart()->legend()->markers(m_endTreatDash).front()->
                setVisible(false);
        if(m_rec){
            m_chartView->chart()->legend()->markers(m_recDash).front()->
                    setVisible(false);
        }
        break;
    }

    case 8:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_cVegfStat->addSeries(m_sDash);
        m_cVegfStat->addSeries(m_endTreatDash);
        m_sDash->attachAxis(m_xVegfStat);
        m_sDash->attachAxis(m_yVegfStat);
        m_endTreatDash->attachAxis(m_xVegfStat);
        m_endTreatDash->attachAxis(m_yVegfStat);
        if(m_rec){
            m_chartView->chart()->removeSeries(m_recDash);
            m_cVegfStat->addSeries(m_recDash);
            m_recDash->attachAxis(m_xVegfStat);
            m_recDash->attachAxis(m_yVegfStat);
        }
        m_chartView->setChart(m_cVegfStat);
        m_chartView->chart()->legend()->markers(m_sDash).front()->
                setVisible(false);
        m_chartView->chart()->legend()->markers(m_endTreatDash).front()->
                setVisible(false);
        if(m_rec){
            m_chartView->chart()->legend()->markers(m_recDash).front()->
                    setVisible(false);
        }
        break;
    }
    }
    m_chartView->setRenderHint(QPainter::Antialiasing);
}


/*------------------------------------------------------------------------------
 * This slot changes the current map iteration.
 * Inputs:
 *  - iter: map iteration to be visualised.
------------------------------------------------------------------------------*/

void OutWindow::changeIter(const int iter){
    drawMap(m_selMap->currentIndex(), iter);
    drawChartDashLine(iter);
}


/*------------------------------------------------------------------------------
 * This slot changes the current map.
 * Inputs:
 *  - numMap: map to be visualised (0, state; 1, cell cycle; 2, pO2; 3, VEGF).
------------------------------------------------------------------------------*/

void OutWindow::changeNumMap(const int numMap){
    drawMap(numMap, m_slider->value());
}


/*------------------------------------------------------------------------------
 * This function draws the chart dash line at the current iteration.
 * Inputs:
 *  - iter: current iteration.
------------------------------------------------------------------------------*/

void OutWindow::drawChartDashLine(const int iter){
    QVector<QPointF> points;
    points.push_back(QPointF(iter, 0.0));
    points.push_back(QPointF(iter, 100.0));
    m_sDash->replace(points);
}


/*------------------------------------------------------------------------------
 * This function draws the current map.
 * Inputs:
 *  - numMap: map to be visualised (0, state; 1, cell cycle; 2, pO2; 3, VEGF),
 *  - mapIter: current iteration.
------------------------------------------------------------------------------*/

void OutWindow::drawMap(const int numMap, const int mapIter){
    const int iter(mapIter / m_simTimeStep);
    const double normPO2(240.0 / m_maxpO2);
    const double normVegf(240.0 / m_maxvegf);
    int colour;
    QColor colourQ;

    switch(numMap){
    case 0:{
        int iSF, jSF;
        for(int i(0); i < m_mapNrow; i++){
            iSF = i * m_mapSclFac;
            for(int j(0); j < m_mapNcol; j++){
                jSF = j * m_mapSclFac;
                colour = m_state[iter][i * m_mapNcol + j] - 1;
                for(int si(0); si < m_mapSclFac; si++){
                    for(int sj(0); sj < m_mapSclFac; sj++){
                        m_imState->setPixel(jSF + sj, iSF + si, colour);
                    }
                }
            }
        }
        m_pixState->convertFromImage(*m_imState);
        m_map->setPixmap(m_pixState->scaledToHeight(500));
        m_legendSt->show();
        m_legendCyc->hide();
        m_legendPO2->hide();
        m_legendVegf->hide();
        break;
    }

    case 1:{
        int iSF, jSF;
        for(int i(0); i < m_mapNrow; i++){
            iSF = i * m_mapSclFac;
            for(int j(0); j < m_mapNcol; j++){
                jSF = j * m_mapSclFac;
                colour = m_timer[iter][i * m_mapNcol + j];
                for(int si(0); si < m_mapSclFac; si++){
                    for(int sj(0); sj < m_mapSclFac; sj++){
                        m_imTimer->setPixel(jSF + sj, iSF + si, colour);
                    }
                }
            }
        }

        m_pixTimer->convertFromImage(*m_imTimer);
        m_map->setPixmap(m_pixTimer->scaledToHeight(500));
        m_legendSt->hide();
        m_legendCyc->show();
        m_legendPO2->hide();
        m_legendVegf->hide();
        break;
    }

    case 2:{
        int iSF, jSF;
        for(int i(0); i < m_mapNrow; i++){
            iSF = i * m_mapSclFac;
            for(int j(0); j < m_mapNcol; j++){
                jSF = j * m_mapSclFac;
                colour = normPO2 * (m_maxpO2 - m_pO2[iter][i * m_mapNcol + j]);
                colourQ = QColor::fromHsv(colour, 200, 255);
                for(int si(0); si < m_mapSclFac; si++){
                    for(int sj(0); sj < m_mapSclFac; sj++){
                        m_imPO2->setPixelColor(jSF + sj, iSF + si, colourQ);
                    }
                }
            }
        }

        m_pixPO2->convertFromImage(*m_imPO2);
        m_map->setPixmap(m_pixPO2->scaledToHeight(500));
        m_legendCyc->hide();
        m_legendSt->hide();
        m_legendPO2->show();
        m_legendVegf->hide();
        break;
    }

    case 3:{
        int iSF, jSF;
        for(int i(0); i < m_mapNrow; i++){
            iSF = i * m_mapSclFac;
            for(int j(0); j < m_mapNcol; j++){
                jSF = j * m_mapSclFac;
                colour = normVegf * (m_maxvegf -
                                     m_vegf[iter][i * m_mapNcol + j]);
                colourQ = QColor::fromHsv(colour, 200, 255);
                for(int si(0); si < m_mapSclFac; si++){
                    for(int sj(0); sj < m_mapSclFac; sj++){
                        m_imVegf->setPixelColor(jSF + sj, iSF + si, colourQ);
                    }
                }
            }
        }

        m_pixVegf->convertFromImage(*m_imVegf);
        m_map->setPixmap(m_pixVegf->scaledToHeight(500));
        m_legendCyc->hide();
        m_legendSt->hide();
        m_legendPO2->hide();
        m_legendVegf->show();
        break;
    }
    }
}


/*------------------------------------------------------------------------------
 * This slot launches a simulation with the same parameters.
 * TO-DO
------------------------------------------------------------------------------*/

void OutWindow::newSim(){
    new InWindow;
    close();
}


/*------------------------------------------------------------------------------
 * This slot plays or pauses the evolution of the current map.
------------------------------------------------------------------------------*/

void OutWindow::play(){
    m_stPlay = !m_stPlay;
    m_play->setIcon(QIcon("../Figures/pause.png"));
    for(int k(m_timeS->value()); k <= m_timeS->maximum(); k++){
        if(!m_stPlay){
            m_play->setIcon(QIcon("../Figures/play.png"));
            return;
        }
        drawMap(m_selMap->currentIndex(), k);
        emit updateSlider(k);
        qApp->processEvents();
    }
    m_stPlay = !m_stPlay;
    m_play->setIcon(QIcon("../Figures/play.png"));
}


/*------------------------------------------------------------------------------
 * This slot saves all the iterations of the current map.
------------------------------------------------------------------------------*/

void OutWindow::saveAllMaps(){
    QString fileName;
    switch (m_selMap->currentIndex()){
    case 0:{
        QDir().mkdir("../Figures/state");
        QDir().mkdir("../Figures/tumDens");
        const int K(m_state.size() * m_simTimeStep);

        for(int k(0); k < K; k += m_simTimeStep){
            drawMap(m_selMap->currentIndex(), k);
            emit updateSlider(k);
            fileName = QString::number(k);
            fileName = fileName.rightJustified(4, '0');
            fileName = "../Figures/state/" + fileName + ".png";
            m_pixState->save(fileName);
            /*fileName = QString::number(k);
            fileName = fileName.rightJustified(4, '0');
            fileName = "../Figures/tumDens/" + fileName + ".png";
            m_chartView->grab().save(fileName);*/
            qApp->processEvents();
        }
        break;
    }

    case 1:{
        QDir().mkdir("../Figures/timer");
        const int K(m_timer.size() * m_simTimeStep);

        for(int k(0); k < K; k += m_simTimeStep){
            drawMap(m_selMap->currentIndex(), k);
            emit updateSlider(k);
            fileName = QString::number(k);
            fileName = fileName.rightJustified(4, '0');
            fileName = "../Figures/timer/" + fileName + ".png";
            m_pixTimer->save(fileName);
            qApp->processEvents();
        }
        break;
    }

    case 2:{
        QDir().mkdir("../Figures/pO2");
        const int K(m_pO2.size() * m_simTimeStep);

        for(int k(0); k < K; k += m_simTimeStep){
            drawMap(m_selMap->currentIndex(), k);
            emit updateSlider(k);
            fileName = QString::number(k);
            fileName = fileName.rightJustified(4, '0');
            fileName = "../Figures/pO2/" + fileName + ".png";
            m_pixPO2->save(fileName);
            qApp->processEvents();
        }
        break;
    }

    case 3:{
        QDir().mkdir("../Figures/vegf");
        const int K(m_vegf.size() * m_simTimeStep);

        for(int k(0); k < K; k += m_simTimeStep){
            drawMap(m_selMap->currentIndex(), k);
            emit updateSlider(k);
            fileName = QString::number(k);
            fileName = fileName.rightJustified(4, '0');
            fileName = "../Figures/vegf/" + fileName + ".png";
            m_pixVegf->save(fileName);
            qApp->processEvents();
        }
        break;
    }
    }
}


/*------------------------------------------------------------------------------
 * This slot saves the current chart.
------------------------------------------------------------------------------*/

void OutWindow::saveChart(){
    QString fileName = QFileDialog::getSaveFileName(this, "Save chart",
                                                    "../Figures",
                                                    "Images(*.png *.gif *.jpg"
                                                    "*.jpeg)");
    if(!fileName.isEmpty()){
        int nSeries(m_chartView->chart()->series().size());
        m_chartView->chart()->series().at(nSeries - 1)->hide();
        m_chartView->chart()->series().at(nSeries - 2)->hide();
        m_chartView->chart()->series().at(nSeries - 3)->hide();
        m_chartView->grab().save(fileName);
        m_chartView->chart()->series().at(nSeries - 1)->show();
        m_chartView->chart()->series().at(nSeries - 2)->show();
        m_chartView->chart()->series().at(nSeries - 3)->show();
        m_chartView->chart()->legend()->markers().at(nSeries - 1)->
                setVisible(false);
        m_chartView->chart()->legend()->markers().at(nSeries - 2)->
                setVisible(false);
        m_chartView->chart()->legend()->markers().at(nSeries - 3)->
                setVisible(false);
    }
}


/*------------------------------------------------------------------------------
 * This slot saves the current iteration of the current map.
------------------------------------------------------------------------------*/

void OutWindow::saveMap(){
    QString fileName = QFileDialog::getSaveFileName(this, "Save map",
                                                    "../Figures",
                                                    "Images(*.png *.gif *.jpg"
                                                    "*.jpeg)");
    if(!fileName.isEmpty()){
        m_map->pixmap()->save(fileName);
    }
}


/*------------------------------------------------------------------------------
 * This slot saves the output files on the specified directory.
------------------------------------------------------------------------------*/

void OutWindow::saveOutFiles(){
    QString dirName = QFileDialog::getSaveFileName(this, "Save ouptut files",
                                                   "..");

    if(!dirName.isEmpty()){
        QDir newDir;
        newDir.mkdir(dirName);

        QDir outDir("../OutputFilesGUI");
        QStringList filesList;
        filesList = outDir.entryList();

        for(int i(0); i < filesList.size(); i++){
            QFile::copy("../OutputFilesGUI/" + filesList[i], dirName + "/" +
                        filesList[i]);
        }
    }
}

/*------------------------------------------------------------------------------
 * This slot goes back to the StartWindow.
------------------------------------------------------------------------------*/

void OutWindow::sel(){
    new StartWindow;
    close();
}
