/*#include <fstream>

#include <QGridLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include "outWindow3D.hpp"

namespace QDV = QtDataVisualization;

OutWindow3D::OutWindow3D(std::string nFOutData) : QWidget(){
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
    m_selMapGroup   = new QGroupBox("Select a map", this);
    m_mapGroup      = new QGroupBox(this);

    m_selChart = new QComboBox(m_selChartGroup);
    m_selMap   = new QComboBox(m_selMapGroup);

    m_chartView = new QChartView;

    double a, b, c, d, e, f;

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
        m_cTumDens->setMaximumWidth(1000);

        QLineSeries *sTumDens = new QLineSeries(m_cTumDens);

        fTumDens >> a >> b;
        while(!fTumDens.eof()){
            sTumDens->append(a, b);
            fTumDens >> a >> b;
        }
        fTumDens.close();

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

        m_sDash = new QLineSeries(m_cTumDens);
        m_sDash->append(0.0, 0.0);
        m_sDash->append(0.0, 100.0);
        m_sDash->setPen(QPen(Qt::DashLine));
        m_cTumDens->addSeries(m_sDash);
        m_sDash->attachAxis(m_xTumDens);
        m_sDash->attachAxis(m_yTumDens);

        m_selChart->addItem("Tumor density");
    }

    std::ifstream fVascDens((nFOutData + "/vascDens.res").c_str());

    if(!fVascDens.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening vascDens.res");
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
        QMessageBox::critical(this, "Error", "Problem while opening killedCells.res");
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
        QMessageBox::critical(this, "Error", "Problem while opening hypDens.res");
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
        QMessageBox::critical(this, "Error", "Problem while opening pO2Stat.res");
    }

    else{
        m_cPO2Stat = new QChart;
        m_xPO2Stat = new QValueAxis;
        m_yPO2Stat = new QValueAxis;
        m_cPO2Stat->addAxis(m_xPO2Stat, Qt::AlignBottom);
        m_cPO2Stat->addAxis(m_yPO2Stat, Qt::AlignLeft);
        m_cPO2Stat->setMaximumWidth(1000);

        QLineSeries *sPO2Med = new QLineSeries(m_cPO2Stat);
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
        m_cPO2Stat->setTitle("Evolution of the PO2 statistics");
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

        m_selChart->addItem("PO2 statistics");
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

    m_map = new QDV::Q3DScatter;

    m_map->activeTheme()->setType(QDV::Q3DTheme::ThemeQt);
    m_map->themes().first()->setGridEnabled(false);
    m_map->themes().first()->setLabelTextColor(m_white);
    m_map->themes().first()->setLabelBorderEnabled(false);
    m_map->setShadowQuality(QDV::QAbstract3DGraph::ShadowQualityNone);
    m_map->setReflection(false);
    m_map->scene()->activeCamera()->setCameraPreset(QDV::Q3DCamera::CameraPresetFront);

    std::ifstream fTissueDim((nFOutData + "/tissueDim.dat").c_str());

    if(!fTissueDim.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening tissueDim.dat");
    }

    else{
        fTissueDim >> m_mapNrow >> m_mapNcol >> m_mapNlayer;
        fTissueDim.close();
    }

    int nrowNcolNlayer(m_mapNrow * m_mapNcol * m_mapNlayer);

    m_mapSclFac = 16;

    std::ifstream fSimParam((nFOutData + "/simParam.dat").c_str());

    if(!fSimParam.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening simParam.dat");
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
            for(int j(0); j < nrowNcolNlayer; j++){
                m_state[i].push_back(tempi);
                fState >> tempi;
            }
            qApp->processEvents();
            i++;
        }
        fState.close();

        m_proxyTum     = new QDV::QScatterDataProxy;
        m_proxyTumDam  = new QDV::QScatterDataProxy;
        m_proxyNormVes = new QDV::QScatterDataProxy;
        m_proxyTumVes  = new QDV::QScatterDataProxy;
        m_proxyHypNec  = new QDV::QScatterDataProxy;
        m_proxyMitCat  = new QDV::QScatterDataProxy;
        m_proxyApop    = new QDV::QScatterDataProxy;

        m_sTum     = new QDV::QScatter3DSeries(m_proxyTum);
        m_sTumDam  = new QDV::QScatter3DSeries(m_proxyTumDam);
        m_sNormVes = new QDV::QScatter3DSeries(m_proxyNormVes);
        m_sTumVes  = new QDV::QScatter3DSeries(m_proxyTumVes);
        m_sHypNec  = new QDV::QScatter3DSeries(m_proxyHypNec);
        m_sMitCat  = new QDV::QScatter3DSeries(m_proxyMitCat);
        m_sApop    = new QDV::QScatter3DSeries(m_proxyApop);

        m_sTum->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
        m_sTum->setMeshSmooth(true);
        m_sTum->setBaseColor(m_blueTum);
        m_sTum->setItemSize(0.1);

        m_sTumDam->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
        m_sTumDam->setMeshSmooth(true);
        m_sTumDam->setBaseColor(m_blueTumDam);
        m_sTumDam->setItemSize(0.1);

        m_sNormVes->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
        m_sNormVes->setMeshSmooth(true);
        m_sNormVes->setBaseColor(m_red);
        m_sNormVes->setItemSize(0.1);

        m_sTumVes->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
        m_sTumVes->setMeshSmooth(true);
        m_sTumVes->setBaseColor(m_redTum);
        m_sTumVes->setItemSize(0.1);

        m_sHypNec->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
        m_sHypNec->setMeshSmooth(true);
        m_sHypNec->setBaseColor(m_yellowNec);
        m_sHypNec->setItemSize(0.1);

        m_sMitCat->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
        m_sMitCat->setMeshSmooth(true);
        m_sMitCat->setBaseColor(m_brownMit);
        m_sMitCat->setItemSize(0.1);

        m_sApop->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
        m_sApop->setMeshSmooth(true);
        m_sApop->setBaseColor(m_greenApop);
        m_sApop->setItemSize(0.1);

        m_datTum     = new QDV::QScatterDataArray;
        m_datTumDam  = new QDV::QScatterDataArray;
        m_datNormVes = new QDV::QScatterDataArray;
        m_datTumVes  = new QDV::QScatterDataArray;
        m_datHypNec  = new QDV::QScatterDataArray;
        m_datMitCat  = new QDV::QScatterDataArray;
        m_datApop    = new QDV::QScatterDataArray;

        m_legendSt = new QGroupBox(m_mapGroup);

        m_fib     = new QLabel("Fibroblasts", m_legendSt);
        m_tum     = new QLabel("Tumor cells", m_legendSt);
        m_tumDam  = new QLabel("Damaged tumor cells", m_legendSt);
        m_normVes = new QLabel("Normal vessels", m_legendSt);
        m_tumVes  = new QLabel("Tumor vessels", m_legendSt);
        m_hypNec  = new QLabel("Hypoxic necrotic cells", m_legendSt);
        m_mitCat  = new QLabel("Mitotic catastrophe cells", m_legendSt);
        m_apop    = new QLabel("Apoptotic cells", m_legendSt);

        QImage fibIm     = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage tumIm     = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage tumDamIm  = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage normVesIm = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage tumVesIm  = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage hypNecIm  = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage mitCatIm  = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage apopIm    = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);

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
            for(int j(0); j < nrowNcolNlayer; j++){
                m_timer[i].push_back(tempi);
                fTimer >> tempi;
            }
            qApp->processEvents();
            i++;
        }
        fTimer.close();

        m_proxyG1 = new QDV::QScatterDataProxy;
        m_proxyS = new QDV::QScatterDataProxy;
        m_proxyG2 = new QDV::QScatterDataProxy;
        m_proxyM = new QDV::QScatterDataProxy;
        m_proxyG0 = new QDV::QScatterDataProxy;

        m_sG1 = new QDV::QScatter3DSeries(m_proxyG1);
        m_sS = new QDV::QScatter3DSeries(m_proxyS);
        m_sG2 = new QDV::QScatter3DSeries(m_proxyG2);
        m_sM = new QDV::QScatter3DSeries(m_proxyM);
        m_sG0 = new QDV::QScatter3DSeries(m_proxyG0);

        m_sG1->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
        m_sG1->setMeshSmooth(true);
        m_sG1->setBaseColor(m_blueG1);
        m_sG1->setItemSize(0.1);

        m_sS->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
        m_sS->setMeshSmooth(true);
        m_sS->setBaseColor(m_green);
        m_sS->setItemSize(0.1);

        m_sG2->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
        m_sG2->setMeshSmooth(true);
        m_sG2->setBaseColor(m_orange);
        m_sG2->setItemSize(0.1);

        m_sM->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
        m_sM->setMeshSmooth(true);
        m_sM->setBaseColor(m_violet);
        m_sM->setItemSize(0.1);

        m_sG0->setItemLabelFormat(QStringLiteral("@xTitle: @xLabel @yTitle: @yLabel @zTitle: @zLabel"));
        m_sG0->setMeshSmooth(true);
        m_sG0->setBaseColor(m_redG0);
        m_sG0->setItemSize(0.1);

        m_datG1 = new QDV::QScatterDataArray;
        m_datS = new QDV::QScatterDataArray;
        m_datG2 = new QDV::QScatterDataArray;
        m_datM = new QDV::QScatterDataArray;
        m_datG0 = new QDV::QScatterDataArray;

        m_legendCyc = new QGroupBox(m_mapGroup);

        m_G1 = new QLabel(QString::fromStdString("G1" + std::string(25, ' ')), m_legendCyc);
        m_S  = new QLabel("S", m_legendCyc);
        m_G2 = new QLabel("G2", m_legendCyc);
        m_M  = new QLabel("M", m_legendCyc);
        m_G0 = new QLabel("G0", m_legendCyc);

        QImage G1Im = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage SIm = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage G2Im = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
        QImage MIm = QImage(m_mapSclFac, m_mapSclFac, QImage::Format_RGB32);
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
    }

    m_time   = new QLabel("Time (h)", m_mapGroup);
    m_timeS  = new QSpinBox(m_mapGroup);
    m_play   = new QPushButton(m_mapGroup);
    m_slider = new QSlider(Qt::Horizontal, this);

    m_timeS->setMinimum(0);
    m_timeS->setMaximum((m_state.size() - 1) * m_simTimeStep);
    m_play->setIcon(QIcon("../Figures/play.png"));
    m_play->setFixedSize(40, 40);
    m_play->setIconSize(QSize(30, 30));
    m_stPlay = false;
    m_slider->setMinimum(0);
    m_slider->setMaximum((m_state.size() - 1) * m_simTimeStep);

    m_sel          = new QPushButton("Select a new mode", this);
    m_change       = new QPushButton("Change input and parameters", this);
    m_saveOutFiles = new QPushButton("Save output files", this);
    m_saveAllMaps  = new QPushButton("Save maps of every iteration", this);
    m_saveChart    = new QPushButton("Save chart", this);
    m_saveMap      = new QPushButton("Save map of current iteration", this);
    m_newSim       = new QPushButton("New simulation", this);

    m_mapCont = new QWidget;
    m_mapCont = QWidget::createWindowContainer(m_map);

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
    mapLayout->addWidget(m_mapCont, 0, 0);
    mapLayout->addWidget(m_legendSt, 0, 1);
    mapLayout->addWidget(m_legendCyc, 0, 1);
    //mapLayout->addWidget(m_legendOxy, 0, 1);
    mapLayout->addWidget(m_slider, 1, 0);
    mapLayout->addLayout(timeLayout, 2, 0);
    mapLayout->addWidget(m_play, 1, 1, 2, 1, Qt::AlignCenter);
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

    setWindowTitle("Radiotherapy Simulator");
    setWindowIcon(QIcon("../Figures/logo.png"));

    QObject::connect(m_selMap, SIGNAL(currentIndexChanged(int)), this, SLOT(changeNumMap(int)));
    QObject::connect(m_selChart, SIGNAL(currentIndexChanged(int)), this, SLOT(changeChart(int)));
    QObject::connect(m_slider, SIGNAL(valueChanged(int)), m_timeS, SLOT(setValue(int)));
    QObject::connect(m_timeS, SIGNAL(valueChanged(int)), m_slider, SLOT(setValue(int)));
    QObject::connect(m_timeS, SIGNAL(valueChanged(int)), this, SLOT(changeIter(int)));
    QObject::connect(this, SIGNAL(updateSlider(int)), m_slider, SLOT(setValue(int)));
    QObject::connect(m_play, SIGNAL(clicked()), this, SLOT(play()));
    QObject::connect(m_sel, SIGNAL(clicked()), this, SLOT(sel()));
    QObject::connect(m_change, SIGNAL(clicked()), this, SLOT(change()));
    QObject::connect(m_newSim, SIGNAL(clicked()),this, SLOT(newSim()));
    QObject::connect(m_saveOutFiles, SIGNAL(clicked()), this, SLOT(saveOutFiles()));
    QObject::connect(m_saveChart, SIGNAL(clicked()), this, SLOT(saveChart()));
    QObject::connect(m_saveMap, SIGNAL(clicked()), this, SLOT(saveMap()));
    QObject::connect(m_saveAllMaps, SIGNAL(clicked()), this, SLOT(saveAllMaps()));

    drawMap(m_selMap->currentIndex(), m_slider->value());
    showMaximized();
    resize(1300, 1000);
    setGeometry(QStyle::alignedRect(Qt::LeftToRight, Qt::AlignCenter, size(),
                                    qApp->desktop()->availableGeometry()));
    show();
}


void OutWindow3D::change(){
    new InWindow;
    close();
}


void OutWindow3D::changeChart(int numChart){
    switch(numChart){
    case 0:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_cTumDens->addSeries(m_sDash);
        m_sDash->attachAxis(m_xTumDens);
        m_sDash->attachAxis(m_yTumDens);
        m_chartView->setChart(m_cTumDens);
        break;
    }
    case 1:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_cVascDens->addSeries(m_sDash);
        m_sDash->attachAxis(m_xVascDens);
        m_sDash->attachAxis(m_yVascDens);
        m_chartView->setChart(m_cVascDens);
        m_chartView->chart()->legend()->markers(m_sDash).front()->setVisible(false);
        break;
    }
    case 2:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_cKilledCells->addSeries(m_sDash);
        m_sDash->attachAxis(m_xKilledCells);
        m_sDash->attachAxis(m_yKilledCells);
        m_chartView->setChart(m_cKilledCells);
        break;
    }
    case 3:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_cCycle->addSeries(m_sDash);
        m_sDash->attachAxis(m_xCycle);
        m_sDash->attachAxis(m_yCycle);
        m_chartView->setChart(m_cCycle);
        m_chartView->chart()->legend()->markers(m_sDash).front()->setVisible(false);
        break;
    }
    case 4:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_cHypDens->addSeries(m_sDash);
        m_sDash->attachAxis(m_xHypDens);
        m_sDash->attachAxis(m_yHypDens);
        m_chartView->setChart(m_cHypDens);
        break;
    }
    case 5:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_cPO2Stat->addSeries(m_sDash);
        m_sDash->attachAxis(m_xPO2Stat);
        m_sDash->attachAxis(m_yPO2Stat);
        m_chartView->setChart(m_cPO2Stat);
        m_chartView->chart()->legend()->markers(m_sDash).front()->setVisible(false);
        break;
    }
    case 6:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_cVegfStat->addSeries(m_sDash);
        m_sDash->attachAxis(m_xVegfStat);
        m_sDash->attachAxis(m_yVegfStat);
        m_chartView->setChart(m_cVegfStat);
        m_chartView->chart()->legend()->markers(m_sDash).front()->setVisible(false);
        break;
    }
    }
    m_chartView->setRenderHint(QPainter::Antialiasing);
}


void OutWindow3D::changeIter(int iter){
    drawMap(m_selMap->currentIndex(), iter);
    drawChartDashLine(iter);
}


void OutWindow3D::changeNumMap(int numMap){
    drawMap(numMap, m_slider->value());
}


void OutWindow3D::drawChartDashLine(int iter){
    QVector<QPointF> points;
    points.push_back(QPointF(iter, 0.0));
    points.push_back(QPointF(iter, 100.0));
    m_sDash->replace(points);
}


void OutWindow3D::drawMap(int numMap, int mapIter){
    while(m_map->seriesList().size()){
        m_map->removeSeries(m_map->seriesList().back());
    }
    int nrowNcol(m_mapNrow * m_mapNcol);

    switch(numMap){
    case 0:{
        m_datTum->clear();
        m_datTumDam->clear();
        m_datNormVes->clear();
        m_datTumVes->clear();
        m_datHypNec->clear();
        m_datMitCat->clear();
        m_datApop->clear();

        for(int i(0); i < m_mapNrow; i++){
            for(int j(0); j < m_mapNcol; j++){
                for(int l(0); l < m_mapNlayer; l++){
                    if(m_state[mapIter / m_simTimeStep][l * nrowNcol + i * m_mapNcol + j] == 2){
                        m_datTum->push_back(QVector3D(i, l, j));
                    }
                    if(m_state[mapIter / m_simTimeStep][l * nrowNcol + i * m_mapNcol + j] == 3){
                        m_datTumDam->push_back(QVector3D(i, l, j));
                    }
                    else if(m_state[mapIter / m_simTimeStep][l * nrowNcol + i * m_mapNcol + j] == 4){
                        m_datNormVes->push_back(QVector3D(i, l, j));
                    }
                    else if(m_state[mapIter / m_simTimeStep][l * nrowNcol + i * m_mapNcol + j] == 5){
                        m_datTumVes->push_back(QVector3D(i, l, j));
                    }
                    else if(m_state[mapIter / m_simTimeStep][l * nrowNcol + i * m_mapNcol + j] == 6){
                        m_datHypNec->push_back(QVector3D(i, l, j));
                    }
                    else if(m_state[mapIter / m_simTimeStep][l * nrowNcol + i * m_mapNcol + j] == 7){
                        m_datMitCat->push_back(QVector3D(i, l, j));
                    }
                    else if(m_state[mapIter / m_simTimeStep][l * nrowNcol + i * m_mapNcol + j] == 8){
                        m_datApop->push_back(QVector3D(i, l, j));
                    }
                }
            }
        }
        m_map->addSeries(m_sTum);
        m_map->addSeries(m_sTumDam);
        m_map->addSeries(m_sNormVes);
        m_map->addSeries(m_sTumVes);
        m_map->addSeries(m_sHypNec);
        m_map->addSeries(m_sMitCat);
        m_map->addSeries(m_sApop);
        m_map->seriesList().at(0)->dataProxy()->resetArray(m_datTum);
        m_map->seriesList().at(1)->dataProxy()->resetArray(m_datTumDam);
        m_map->seriesList().at(2)->dataProxy()->resetArray(m_datNormVes);
        m_map->seriesList().at(3)->dataProxy()->resetArray(m_datTumVes);
        m_map->seriesList().at(4)->dataProxy()->resetArray(m_datHypNec);
        m_map->seriesList().at(5)->dataProxy()->resetArray(m_datMitCat);
        m_map->seriesList().at(6)->dataProxy()->resetArray(m_datApop);
        m_legendSt->show();
        m_legendCyc->hide();
        //m_legendOxy->hide();
        break;
    }

    case 1:{
        m_datG1->clear();
        m_datS->clear();
        m_datG2->clear();
        m_datM->clear();
        m_datG0->clear();

        for(int i(0); i < m_mapNrow; i++){
            for(int j(0); j < m_mapNcol; j++){
                for(int l(0); l < m_mapNlayer; l++){
                    if(m_timer[mapIter / m_simTimeStep][l * nrowNcol + i * m_mapNcol + j] == 1){
                        m_datG1->push_back(QVector3D(i, l, j));
                    }
                    if(m_timer[mapIter / m_simTimeStep][l * nrowNcol + i * m_mapNcol + j] == 2){
                        m_datS->push_back(QVector3D(i, l, j));
                    }
                    else if(m_timer[mapIter / m_simTimeStep][l * nrowNcol + i * m_mapNcol + j] == 3){
                        m_datG2->push_back(QVector3D(i, l, j));
                    }
                    else if(m_timer[mapIter / m_simTimeStep][l * nrowNcol + i * m_mapNcol + j] == 4){
                        m_datM->push_back(QVector3D(i, l, j));
                    }
                    else if(m_timer[mapIter / m_simTimeStep][l * nrowNcol + i * m_mapNcol + j] == 5){
                        m_datG0->push_back(QVector3D(i, l, j));
                    }
                }
            }
        }

        m_map->addSeries(m_sG1);
        m_map->addSeries(m_sS);
        m_map->addSeries(m_sG2);
        m_map->addSeries(m_sM);
        m_map->addSeries(m_sG0);
        m_map->seriesList().at(0)->dataProxy()->resetArray(m_datG1);
        m_map->seriesList().at(1)->dataProxy()->resetArray(m_datS);
        m_map->seriesList().at(2)->dataProxy()->resetArray(m_datG2);
        m_map->seriesList().at(3)->dataProxy()->resetArray(m_datM);
        m_map->seriesList().at(4)->dataProxy()->resetArray(m_datG0);
        m_legendCyc->show();
        m_legendSt->hide();
        //m_legendOxy->hide();
        break;
    }

    case 2:{
        m_legendCyc->hide();
        m_legendSt->hide();
        m_legendOxy->show();
        break;
    }
    }
}


void OutWindow3D::newSim(){
    new InWindow;
    close();
}


void OutWindow3D::play(){
    m_stPlay = !m_stPlay;
    m_play->setIcon(QIcon("../Figures/pause.png"));
    for(int k(m_timeS->value()); k < m_timeS->maximum(); k++){
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


void OutWindow3D::saveAllMaps(){
    QString fileName;
    switch (m_selMap->currentIndex()){
    case 0:{
        QDir().mkdir("../Figures/state");

        for(int k(0); k < m_state.size(); k++){
            drawMap(m_selMap->currentIndex(), k);
            emit updateSlider(k);
            fileName = "../Figures/state/" + QString::number(k) + ".png";
            m_mapCont->grab().save(fileName);
            qApp->processEvents();
        }
        break;
    }

    case 1:{
        QDir().mkdir("../Figures/timer");

        for(int k(0); k < m_timer.size(); k++){
            drawMap(m_selMap->currentIndex(), k);
            emit updateSlider(k);
            fileName = "../Figures/timer/" + QString::number(k) + ".png";
            m_mapCont->grab().save(fileName);
            qApp->processEvents();
        }
        break;
    }

    case 2:{
        QDir().mkdir("../Figures/pO2");

        for(int k(0); k < m_pO2.size(); k++){
            drawMap(m_selMap->currentIndex(), k);
            emit updateSlider(k);
            fileName = "../Figures/pO2/" + QString::number(k) + ".png";
            m_mapCont->grab().save(fileName);
            qApp->processEvents();
        }
        break;
    }
    }
}


void OutWindow3D::saveChart(){
    QString fileName = QFileDialog::getSaveFileName(this, "Save chart", "../Figures",
                                                    "Images(*.png *.gif *.jpg *.jpeg)");
    if(!fileName.isEmpty()){
        m_chartView->chart()->series().back()->hide();
        m_chartView->grab().save(fileName);
        m_chartView->chart()->series().back()->show();
        m_chartView->chart()->legend()->markers().back()->setVisible(false);
    }
}


void OutWindow3D::saveMap(){
    QString fileName = QFileDialog::getSaveFileName(this, "Save map", "../Figures",
                                                    "Images(*.png *.gif *.jpg *.jpeg)");
    if(!fileName.isEmpty()){
        m_mapCont->grab().save(fileName);
    }
}


void OutWindow3D::saveOutFiles(){
    QString dirName = QFileDialog::getSaveFileName(this, "Save ouptut files", "..");

    if(!dirName.isEmpty()){
        QDir newDir;
        newDir.mkdir(dirName);

        QDir outDir("../OutputFiles");
        QStringList filesList;
        filesList = outDir.entryList();

        for(int i(0); i < filesList.size(); i++){
            QFile::copy("../OutputFiles/" + filesList[i], dirName + "/" + filesList[i]);
        }
    }
}


void OutWindow3D::sel(){
    new StartWindow;
    close();
}*/
