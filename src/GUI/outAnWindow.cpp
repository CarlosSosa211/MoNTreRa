#include <fstream>
#include <iomanip>
#include <iostream>

#include <QGridLayout>
#include <QHBoxLayout>
#include <QMessageBox>
#include <QVBoxLayout>

#include "outAnWindow.hpp"

OutAnWindow::OutAnWindow(std::string nFOutData) : QWidget(){
    m_white      = QColor(255, 255, 255);
    m_blueTum    = QColor(46, 165, 225);
    m_red        = QColor(192, 0, 0);
    m_redTum     = QColor(153, 0, 51);
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
    m_selChart = new QComboBox(m_selChartGroup);
    m_chartView = new QChartView;

    std::ifstream fSimParam((nFOutData + "/simParam.dat").c_str());

    if(!fSimParam.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening simParam.dat");
    }

    else{
        fSimParam >> m_simTime >> m_simTime >> m_simTimeStep >> m_simNum;;
        fSimParam.close();
    }

    std::ifstream fEndTreatTime((nFOutData + "/endTreatTumDens.res").c_str());

    if(!fEndTreatTime.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening endTreatTumDens.res");
    }

    else{
        fEndTreatTime >> m_endTreatTime;
        fEndTreatTime.close();
    }

    std::ifstream fRecTime((nFOutData + "/rec.res").c_str());

    if(!fRecTime.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening rec.res");
    }

    else{
        fRecTime >> m_recTime;
        fRecTime.close();
    }

    m_endTreatDash = new QLineSeries;
    m_endTreatDash->append(m_endTreatTime, 0.0);
    m_endTreatDash->append(m_endTreatTime, 100.0);
    QPen pen;
    pen.setStyle(Qt::DashLine);
    pen.setColor(m_green);
    m_endTreatDash->setPen(pen);

    m_recDash = new QLineSeries;
    m_recDash->append(m_recTime, 0.0);
    m_recDash->append(m_recTime, 100.0);
    pen.setColor(m_red);
    m_recDash->setPen(pen);

    double a, b[m_simNum], c[m_simNum], d[m_simNum], e[m_simNum], f[m_simNum];

    std::ifstream fTumDens((nFOutData + "/tumDens.res").c_str());

    if(!fTumDens.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening tumDens.res");
    }

    else{
        m_cTumDens = new QChart;
        m_xTumDens = new QValueAxis;
        m_yTumDens = new QValueAxis;
        m_cTumDens->addAxis(m_xTumDens, Qt::AlignBottom);
        m_cTumDens->addAxis(m_yTumDens, Qt::AlignLeft);

        std::vector<QLineSeries *> sTumDens(m_simNum);

        fTumDens >> a;
        for(int j(0); j < m_simNum; j++){
            fTumDens >> b[j];
        }
        while(!fTumDens.eof()){
            for(int j(0); j < m_simNum; j++){
                sTumDens.at(j)->append(a, b[j]);
            }
            fTumDens >> a;
            for(int j(0); j < m_simNum; j++){
                fTumDens >> b[j];
            }
        }
        fTumDens.close();

        for(int j(0); j < m_simNum; j++){
            m_cTumDens->addSeries(sTumDens.at(j));
            sTumDens.at(j)->attachAxis(m_xTumDens);
            sTumDens.at(j)->attachAxis(m_yTumDens);
        }
        m_cTumDens->legend()->hide();
        m_cTumDens->setTitle("Evolution of the tumor density");
        m_xTumDens->setTitleText("Time (h)");
        m_yTumDens->setTitleText("Tumor density");
        m_xTumDens->setLabelFormat("%i");
        m_yTumDens->setLabelFormat("%i");
        m_xTumDens->applyNiceNumbers();
        m_yTumDens->applyNiceNumbers();
        m_yTumDens->setMin(0.0);

        m_cTumDens->addSeries(m_endTreatDash);
        m_endTreatDash->attachAxis(m_xTumDens);
        m_endTreatDash->attachAxis(m_yTumDens);

        m_cTumDens->addSeries(m_recDash);
        m_recDash->attachAxis(m_xTumDens);
        m_recDash->attachAxis(m_yTumDens);

        m_selChart->addItem("Tumor density");
    }

    std::ifstream fTumVol((nFOutData + "/tumVol.res").c_str());

    if(!fTumVol.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening tumVol.res");
    }

    else{
        m_cTumVol = new QChart;
        m_xTumVol = new QValueAxis;
        m_yTumVol = new QValueAxis;
        m_cTumVol->addAxis(m_xTumVol, Qt::AlignBottom);
        m_cTumVol->addAxis(m_yTumVol, Qt::AlignLeft);

        std::vector<QLineSeries *> sTumVol(m_simNum);

        for(int j(0); j < m_simNum; j++){
            fTumVol >> a >> b[j];
        }
        while(!fTumVol.eof()){
            for(int j(0); j < m_simNum; j++){
                sTumVol.at(j)->append(a, b[j]);
            }
            fTumVol >> a;
            for(int j(0); j < m_simNum; j++){
                fTumVol >> b[j];
            }
        }
        fTumVol.close();

        for(int j(0); j < m_simNum; j++){
            m_cTumVol->addSeries(sTumVol.at(j));
            sTumVol.at(j)->attachAxis(m_xTumVol);
            sTumVol.at(j)->attachAxis(m_yTumVol);
        }
        m_cTumVol->legend()->hide();
        m_cTumVol->setTitle("Evolution of the tumor volume");
        m_xTumVol->setTitleText("Time (h)");
        m_yTumVol->setTitleText("Tumor volume (mm³)");
        m_xTumVol->setLabelFormat("%i");
        m_yTumVol->setLabelFormat("%4.2f");
        m_xTumVol->applyNiceNumbers();
        m_yTumVol->applyNiceNumbers();
        m_yTumVol->setMin(0.0);

        m_selChart->addItem("Tumor volume");
    }

    std::ifstream fVascDens((nFOutData + "/vascDens.res").c_str());

    /*if(!fVascDens.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening vascDens.res");
    }

    else{
        m_cVascDens = new QChart;
        m_xVascDens = new QValueAxis;
        m_yVascDens = new QValueAxis;
        m_cVascDens->addAxis(m_xVascDens, Qt::AlignBottom);
        m_cVascDens->addAxis(m_yVascDens, Qt::AlignLeft);

        std::vector <QLineSeries *> sVascDens(m_simNum);
        std::vector <QLineSeries *> sNormVascDens(m_simNum);
        std::vector <QLineSeries *> sTumVascDens(m_simNum);

        fVascDens >> a;
        for(int j(0); j < m_simNum; j++){
            fVascDens >> b[j] >> c[j] >> d[j];
        }
        while(!fVascDens.eof()){
            for(int j(0); j < m_simNum; j++){
                sVascDens.at(j)->append(a, b[j]);
                sNormVascDens.at(j)->append(a, c[j]);
                sTumVascDens.at(j)->append(a, d[j]);
            }
            fVascDens >> a;
            for(int j(0); j < m_simNum; j++){
                fVascDens >> b[j] >> c[j] >> d[j];
            }
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
        QMessageBox::critical(this, "Error", "Problem while opening killedCells.res");
    }

    else{
        m_cKilledCells = new QChart;
        m_xKilledCells = new QValueAxis;
        m_yKilledCells = new QValueAxis;
        m_cKilledCells->addAxis(m_xKilledCells, Qt::AlignBottom);
        m_cKilledCells->addAxis(m_yKilledCells, Qt::AlignLeft);

        QLineSeries *sKilledCells = new QLineSeries(m_cKilledCells);

        fKilledCells >> a >> b;
        while(!fKilledCells.eof()){
            sKilledCells->append(a, b);
            fKilledCells >> a >> b;
        }
        fKilledCells.close();

        m_cKilledCells->addSeries(sKilledCells);
        m_cKilledCells->legend()->hide();
        m_cKilledCells->setTitle("Evolution of the percentage of killed tumor cells");

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
        QMessageBox::critical(this, "Error", "Problem while opening hypDens.res");
    }

    else{
        m_cHypDens = new QChart;
        m_xHypDens = new QValueAxis;
        m_yHypDens = new QValueAxis;
        m_cHypDens->addAxis(m_xHypDens, Qt::AlignBottom);
        m_cHypDens->addAxis(m_yHypDens, Qt::AlignLeft);

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
        QMessageBox::critical(this, "Error", "Problem while opening pO2Stat.res");
    }

    else{
        m_cPO2Stat = new QChart;
        m_xPO2Stat = new QValueAxis;
        m_yPO2Stat = new QValueAxis;
        m_cPO2Stat->addAxis(m_xPO2Stat, Qt::AlignBottom);
        m_cPO2Stat->addAxis(m_yPO2Stat, Qt::AlignLeft);

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
        QMessageBox::critical(this, "Error", "Problem while opening vegfStat.res");
    }

    else{
        m_cVegfStat = new QChart;
        m_xVegfStat = new QValueAxis;
        m_yVegfStat = new QValueAxis;
        m_cVegfStat->addAxis(m_xVegfStat, Qt::AlignBottom);
        m_cVegfStat->addAxis(m_yVegfStat, Qt::AlignLeft);;

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
    }*/

    m_chartView->setChart(m_cTumDens);

    m_sel          = new QPushButton("Select a new mode", this);
    m_change       = new QPushButton("Change input and parameters", this);
    m_saveOutFiles = new QPushButton("Save output files", this);
    m_saveChart    = new QPushButton("Save chart", this);
    m_newSim       = new QPushButton("New simulation", this);

    QVBoxLayout *chartLayout = new QVBoxLayout;
    chartLayout->addWidget(m_selChart);
    m_selChartGroup->setLayout(chartLayout);

    QHBoxLayout *hLayout = new QHBoxLayout;
    hLayout->addWidget(m_sel);
    hLayout->addWidget(m_change);
    hLayout->addWidget(m_saveOutFiles);
    hLayout->addWidget(m_saveChart);
    hLayout->addWidget(m_newSim);

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(m_selChartGroup);
    layout->addWidget(m_chartView);
    layout->addLayout(hLayout);
    setLayout(layout);

    QObject::connect(m_selChart, SIGNAL(currentIndexChanged(int)), this, SLOT(changeChart(int)));
    QObject::connect(m_sel, SIGNAL(clicked()), this, SLOT(sel()));
    QObject::connect(m_change, SIGNAL(clicked()), this, SLOT(change()));
    QObject::connect(m_newSim, SIGNAL(clicked()),this, SLOT(newSim()));
    QObject::connect(m_saveOutFiles, SIGNAL(clicked()), this, SLOT(saveOutFiles()));
    QObject::connect(m_saveChart, SIGNAL(clicked()), this, SLOT(saveChart()));

    setWindowTitle("Radiotherapy Simulator");
    setWindowIcon(QIcon("../Figures/logo.png"));
    showMaximized();
}


void OutAnWindow::change(){
    new InWindow;
    close();
}


void OutAnWindow::changeChart(const int numChart){
    switch(numChart){
    case 0:{
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_chartView->chart()->removeSeries(m_recDash);
        m_cTumDens->addSeries(m_endTreatDash);
        m_cTumDens->addSeries(m_recDash);
        m_endTreatDash->attachAxis(m_xTumDens);
        m_endTreatDash->attachAxis(m_yTumDens);
        m_recDash->attachAxis(m_xTumDens);
        m_recDash->attachAxis(m_yTumDens);
        m_chartView->setChart(m_cTumDens);
        break;
    }

    case 1:{
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_chartView->chart()->removeSeries(m_recDash);
        m_cTumVol->addSeries(m_endTreatDash);
        m_cTumVol->addSeries(m_recDash);
        m_endTreatDash->attachAxis(m_xTumVol);
        m_endTreatDash->attachAxis(m_yTumVol);
        m_recDash->attachAxis(m_xTumVol);
        m_recDash->attachAxis(m_yTumVol);
        m_chartView->setChart(m_cTumVol);
        break;
    }

    case 2:{
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_chartView->chart()->removeSeries(m_recDash);
        m_cVascDens->addSeries(m_endTreatDash);
        m_cVascDens->addSeries(m_recDash);
        m_endTreatDash->attachAxis(m_xVascDens);
        m_endTreatDash->attachAxis(m_yVascDens);
        m_recDash->attachAxis(m_xVascDens);
        m_recDash->attachAxis(m_yVascDens);
        m_chartView->setChart(m_cVascDens);
        m_chartView->chart()->legend()->markers(m_endTreatDash).front()->setVisible(false);
        m_chartView->chart()->legend()->markers(m_recDash).front()->setVisible(false);
        break;
    }

    case 3:{
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_chartView->chart()->removeSeries(m_recDash);
        m_cKilledCells->addSeries(m_endTreatDash);
        m_cKilledCells->addSeries(m_recDash);
        m_endTreatDash->attachAxis(m_xKilledCells);
        m_endTreatDash->attachAxis(m_yKilledCells);
        m_recDash->attachAxis(m_xKilledCells);
        m_recDash->attachAxis(m_yKilledCells);
        m_chartView->setChart(m_cKilledCells);
        break;
    }

    case 4:{
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_chartView->chart()->removeSeries(m_recDash);
        m_cCycle->addSeries(m_endTreatDash);
        m_cCycle->addSeries(m_recDash);
        m_endTreatDash->attachAxis(m_xCycle);
        m_endTreatDash->attachAxis(m_yCycle);
        m_recDash->attachAxis(m_xCycle);
        m_recDash->attachAxis(m_yCycle);
        m_chartView->setChart(m_cCycle);
        m_chartView->chart()->legend()->markers(m_endTreatDash).front()->setVisible(false);
        m_chartView->chart()->legend()->markers(m_recDash).front()->setVisible(false);
        break;
    }

    case 5:{
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_chartView->chart()->removeSeries(m_recDash);
        m_cHypDens->addSeries(m_endTreatDash);
        m_cHypDens->addSeries(m_recDash);
        m_endTreatDash->attachAxis(m_xHypDens);
        m_endTreatDash->attachAxis(m_yHypDens);
        m_recDash->attachAxis(m_xHypDens);
        m_recDash->attachAxis(m_yHypDens);
        m_chartView->setChart(m_cHypDens);
        break;
    }

    case 6:{
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_chartView->chart()->removeSeries(m_recDash);
        m_cPO2Stat->addSeries(m_endTreatDash);
        m_cPO2Stat->addSeries(m_recDash);
        m_endTreatDash->attachAxis(m_xPO2Stat);
        m_endTreatDash->attachAxis(m_yPO2Stat);
        m_recDash->attachAxis(m_xPO2Stat);
        m_recDash->attachAxis(m_yPO2Stat);
        m_chartView->setChart(m_cPO2Stat);
        m_chartView->chart()->legend()->markers(m_endTreatDash).front()->setVisible(false);
        m_chartView->chart()->legend()->markers(m_recDash).front()->setVisible(false);
        break;
    }

    case 7:{
        m_chartView->chart()->removeSeries(m_endTreatDash);
        m_chartView->chart()->removeSeries(m_recDash);
        m_cVegfStat->addSeries(m_endTreatDash);
        m_cVegfStat->addSeries(m_recDash);
        m_endTreatDash->attachAxis(m_xVegfStat);
        m_endTreatDash->attachAxis(m_yVegfStat);
        m_recDash->attachAxis(m_xVegfStat);
        m_recDash->attachAxis(m_yVegfStat);
        m_chartView->setChart(m_cVegfStat);
        m_chartView->chart()->legend()->markers(m_endTreatDash).front()->setVisible(false);
        m_chartView->chart()->legend()->markers(m_recDash).front()->setVisible(false);
        break;
    }
    }
    m_chartView->setRenderHint(QPainter::Antialiasing);
}


void OutAnWindow::newSim(){
    new InAnWindow;
    close();
}


void OutAnWindow::saveChart(){
    QString fileName = QFileDialog::getSaveFileName(this, "Savem_chart", "../Figures",
                                                    "Images(*.png *.gif *.jpg *.jpeg)");
    if(!fileName.isEmpty()){
        int nSeries(m_chartView->chart()->series().size());
        m_chartView->chart()->series().at(nSeries - 1)->hide();
        m_chartView->chart()->series().at(nSeries - 2)->hide();
        m_chartView->grab().save(fileName);
        m_chartView->chart()->series().at(nSeries - 1)->show();
        m_chartView->chart()->series().at(nSeries - 2)->show();
        m_chartView->chart()->series().at(nSeries - 3)->show();
        m_chartView->chart()->legend()->markers().at(nSeries - 1)->setVisible(false);
        m_chartView->chart()->legend()->markers().at(nSeries - 2)->setVisible(false);
    }
}



void OutAnWindow::saveOutFiles(){
    QString dirName = QFileDialog::getSaveFileName(this, "Save ouptut files", "..");

    if(!dirName.isEmpty()){
        QDir newDir;
        newDir.mkdir(dirName);

        QDir outDir("../OutputFilesGUI");
        QStringList filesList;
        filesList = outDir.entryList();

        for(int i(0); i < filesList.size(); i++){
            QFile::copy("../OutputFilesGUI/" + filesList[i], dirName + "/" + filesList[i]);
        }
    }
}

void OutAnWindow::sel(){
    new StartWindow;
    close();
}
