#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <time.h>

#include <QHBoxLayout>
#include <QVBoxLayout>

#include "inAnWindow.hpp"

InAnWindow::InAnWindow() : QWidget(){
    m_param = new QGroupBox("Select the input or parameter to be analyzed and enter the corresponding values",
                            this);
    m_tumDens   = new QCheckBox("Initial tumor density", m_param);
    m_sigmaTum  = new QCheckBox("Standard deviation\nof tumor cells", m_param);
    m_vascDens  = new QCheckBox("Initial vascular density", m_param);
    m_sigmaVasc = new QCheckBox("Standard deviation\nof vessels", m_param);

    m_tumGrowth = new QCheckBox("Tumor growth", m_param);
    m_doubTime  = new QCheckBox("Doubling time of tumor cells (h)", m_param);
    m_edgeOrder = new QCheckBox("Edge order", m_param);

    m_res         = new QCheckBox("Fibroblasts proliferation", m_param);
    m_fibDoubTime = new QCheckBox("Doubling time of fibroblasts (h)", m_param);

    m_ang       = new QCheckBox("Angiogenesis", m_param);
    m_angTime   = new QCheckBox("Doubling time of vessels (h)", m_param);
    m_Dvegf     = new QCheckBox("Diffusion coefficient (μm²/ms)", m_param);
    m_VmaxVegf  = new QCheckBox("Vmax (mol/μm³ms)", m_param);
    m_KmVegf    = new QCheckBox("K (mol/μm³)", m_param);
    m_hypVegf   = new QCheckBox("Hypoxic cells VEGF (mol/μm³)", m_param);
    m_vegfThres = new QCheckBox("VEGF threshold (mol/μm³)", m_param);

    m_treat      = new QCheckBox("Treatment", m_param);
    m_phaseRS    = new QCheckBox("Phase-dependent radiosensitivity", m_param);
    m_alphaGroup = new QGroupBox("Alpha (Gy⁻¹)", m_param);
    m_betaGroup  = new QGroupBox("Beta (Gy⁻²)", m_param);
    m_alpha.push_back(new QCheckBox("Fibroblasts", m_param));
    m_alpha.push_back(new QCheckBox("Tumor cells in phase G1", m_alphaGroup));
    m_alpha.push_back(new QCheckBox("Tumor cells in phase S", m_alphaGroup));
    m_alpha.push_back(new QCheckBox("Tumor cells in phase G2", m_alphaGroup));
    m_alpha.push_back(new QCheckBox("Tumor cells in phase M", m_alphaGroup));
    m_alpha.push_back(new QCheckBox("Tumor cells in phase G0", m_alphaGroup));
    m_alpha.push_back(new QCheckBox("Pre-existing vessels", m_alphaGroup));
    m_alpha.push_back(new QCheckBox("Neocreated vessels", m_alphaGroup));
    m_beta.push_back(new QCheckBox("Fibroblasts", m_betaGroup));
    m_beta.push_back(new QCheckBox("Tumor cells in phase G1", m_betaGroup));
    m_beta.push_back(new QCheckBox("Tumor cells in phase S", m_betaGroup));
    m_beta.push_back(new QCheckBox("Tumor cells in phase G2", m_betaGroup));
    m_beta.push_back(new QCheckBox("Tumor cells in phase M", m_betaGroup));
    m_beta.push_back(new QCheckBox("Tumor cells in phase G0", m_betaGroup));
    m_beta.push_back(new QCheckBox("Pre-existing vessels", m_betaGroup));
    m_beta.push_back(new QCheckBox("Neocreated vessels", m_betaGroup));

    m_arrestTime = new QCheckBox("Time of cycle arrest (h)", m_param);
    m_doseThres  = new QCheckBox("Mitotic catastrophe \ndose threshold (Gy)", m_param);
    m_fraction   = new QCheckBox("Fraction (Gy)", m_param);

    m_oxy         = new QCheckBox("Oxygenation", m_param);
    m_D           = new QCheckBox("Diffusion coefficient (μm²/ms)", m_param);
    m_Vmax        = new QCheckBox("Vmax (mmHg/ms)", m_param);
    m_Km          = new QCheckBox("K (mmHg)", m_param);
    m_pO2NormVes  = new QCheckBox("Normal vessels pO2 (mmHg)", m_param);
    m_pO2TumVes   = new QCheckBox("Tumor vessels pO2 (mmHg)", m_param);
    m_hypThres    = new QCheckBox("Hypoxia threshold (mmHg)", m_param);
    m_hypNecThres = new QCheckBox("Hypoxic necrosis threshold (mmHg)", m_param);

    m_tumDensS   = new QDoubleSpinBox(m_param);
    m_sigmaTumS  = new QDoubleSpinBox(m_param);
    m_vascDensS  = new QDoubleSpinBox(m_param);
    m_sigmaVascS = new QDoubleSpinBox(m_param);

    m_tumGrowthS = new Switch(m_param);
    m_doubTimeS  = new QDoubleSpinBox(m_param);
    m_edgeOrderS = new QSpinBox(m_param);

    m_resS         = new Switch(m_param);
    m_fibDoubTimeS = new QDoubleSpinBox( m_param);

    m_angS       = new Switch(m_param);
    m_angTimeS   = new QDoubleSpinBox(m_param);
    m_DvegfS     = new QDoubleSpinBox(m_param);
    m_VmaxVegfS  = new QDoubleSpinBox(m_param);
    m_KmVegfS    = new QDoubleSpinBox(m_param);
    m_hypVegfS   = new QDoubleSpinBox(m_param);
    m_vegfThresS = new QDoubleSpinBox(m_param);

    m_treatS   = new Switch(m_param);
    m_phaseRSS = new Switch(m_param);
    m_alphaS.push_back(new QDoubleSpinBox(m_param));
    m_alphaS.push_back(new QDoubleSpinBox(m_alphaGroup));
    m_alphaS.push_back(new QDoubleSpinBox(m_alphaGroup));
    m_alphaS.push_back(new QDoubleSpinBox(m_alphaGroup));
    m_alphaS.push_back(new QDoubleSpinBox(m_alphaGroup));
    m_alphaS.push_back(new QDoubleSpinBox(m_alphaGroup));
    m_alphaS.push_back(new QDoubleSpinBox(m_alphaGroup));
    m_alphaS.push_back(new QDoubleSpinBox(m_alphaGroup));
    m_betaS.push_back(new QDoubleSpinBox(m_betaGroup));
    m_betaS.push_back(new QDoubleSpinBox(m_betaGroup));
    m_betaS.push_back(new QDoubleSpinBox(m_betaGroup));
    m_betaS.push_back(new QDoubleSpinBox(m_betaGroup));
    m_betaS.push_back(new QDoubleSpinBox(m_betaGroup));
    m_betaS.push_back(new QDoubleSpinBox(m_betaGroup));
    m_betaS.push_back(new QDoubleSpinBox(m_betaGroup));
    m_betaS.push_back(new QDoubleSpinBox(m_betaGroup));

    m_arrestTimeS = new QDoubleSpinBox(m_param);
    m_doseThresS  = new QDoubleSpinBox(m_param);
    m_fractionS   = new QDoubleSpinBox(m_param);

    m_oxyS         = new Switch(m_param);
    m_DS           = new QDoubleSpinBox(m_param);
    m_VmaxS        = new QDoubleSpinBox(m_param);
    m_KmS          = new QDoubleSpinBox(m_param);
    m_pO2NormVesS  = new QDoubleSpinBox(m_param);
    m_pO2TumVesS   = new QDoubleSpinBox(m_param);
    m_hypThresS    = new QDoubleSpinBox(m_param);
    m_hypNecThresS = new QDoubleSpinBox(m_param);

    m_back       = new QPushButton("Back", this);
    m_cancel     = new QPushButton("Cancel", this);
    m_loadInData = new QPushButton("Load input data", this);
    m_start      = new QPushButton("Start analysis", this);
    m_progress   = new QProgressBar(this);

    m_tumDensS->setMaximum(1.0);
    m_tumDensS->setSingleStep(0.01);
    m_sigmaTumS->setMaximum(2.0);
    m_sigmaTumS->setSingleStep(0.01);
    m_vascDensS->setMaximum(1.0);
    m_vascDensS->setSingleStep(0.01);
    m_sigmaVascS->setMaximum(2.0);
    m_sigmaVascS->setSingleStep(0.01);

    m_tumGrowthS->setChecked(true);
    m_doubTimeS->setMaximum(9999.9);
    m_doubTimeS->setSingleStep(1.0);
    m_doubTimeS->setDecimals(1);
    m_edgeOrderS->setMaximum(10);

    m_resS->setChecked(true);
    m_fibDoubTimeS->setMaximum(9999.9);
    m_fibDoubTimeS->setSingleStep(1.0);
    m_fibDoubTimeS->setDecimals(1);

    m_angS->setChecked(true);
    m_angTimeS->setMaximum(9999.9);
    m_angTimeS->setSingleStep(1.0);
    m_angTimeS->setDecimals(1);
    m_DvegfS->setMaximum(9.9);
    m_DvegfS->setSingleStep(0.1);
    m_DvegfS->setDecimals(1);
    m_VmaxVegfS->setDecimals(3);
    m_KmVegfS->setDecimals(3);
    m_hypVegfS->setDecimals(1);
    m_vegfThresS->setDecimals(1);

    m_treatS->setChecked(true);
    m_phaseRSS->setChecked(true);
    for(int i(0); i < 8; i++){
        m_alphaS.at(i)->setSingleStep(0.001);
        m_alphaS.at(i)->setDecimals(3);
        m_betaS.at(i)->setSingleStep(0.001);
        m_betaS.at(i)->setDecimals(4);
    }
    m_fractionS->setDecimals(1);

    m_oxyS->setChecked(true);
    m_DS->setMaximum(9.9);
    m_DS->setSingleStep(0.1);
    m_DS->setDecimals(1);
    m_VmaxS->setDecimals(3);
    m_KmS->setDecimals(3);
    m_hypNecThresS->setDecimals(3);

    m_progress->setMinimum(0);


    QGridLayout *alphaLayout = new QGridLayout;
    for(int i(0); i < 8; i++){
        alphaLayout->addWidget(m_alpha.at(i), i, 0);
        alphaLayout->addWidget(m_alphaS.at(i), i, 1);
    }
    m_alphaGroup->setLayout(alphaLayout);

    QGridLayout *betaLayout = new QGridLayout;
    for(int i(0); i < 8; i++){
        betaLayout->addWidget(m_beta.at(i), i, 0);
        betaLayout->addWidget(m_betaS.at(i), i, 1);
    }
    m_betaGroup->setLayout(betaLayout);

    QGridLayout *paramLLayout = new QGridLayout;
    paramLLayout->addWidget(m_tumDens, 1, 0);
    paramLLayout->addWidget(m_sigmaTum, 2, 0);
    paramLLayout->addWidget(m_vascDens, 3, 0);
    paramLLayout->addWidget(m_sigmaVasc, 4, 0);
    paramLLayout->addWidget(m_tumGrowth, 5, 0);
    paramLLayout->addWidget(m_doubTime, 6, 0);
    paramLLayout->addWidget(m_edgeOrder, 7, 0);
    paramLLayout->addWidget(m_res, 8, 0);
    paramLLayout->addWidget(m_fibDoubTime, 9, 0);
    paramLLayout->addWidget(m_ang, 10, 0);
    paramLLayout->addWidget(m_angTime, 11, 0);
    paramLLayout->addWidget(m_Dvegf, 12, 0);
    paramLLayout->addWidget(m_VmaxVegf, 13, 0);
    paramLLayout->addWidget(m_KmVegf, 14, 0);
    paramLLayout->addWidget(m_hypVegf, 15, 0);
    paramLLayout->addWidget(m_vegfThres, 16, 0);

    paramLLayout->addWidget(m_tumDensS, 1, 1);
    paramLLayout->addWidget(m_sigmaTumS, 2, 1);
    paramLLayout->addWidget(m_vascDensS, 3, 1);
    paramLLayout->addWidget(m_sigmaVascS, 4, 1);
    paramLLayout->addWidget(m_tumGrowthS, 5, 1);
    paramLLayout->addWidget(m_doubTimeS, 6, 1);
    paramLLayout->addWidget(m_edgeOrderS, 7, 1);
    paramLLayout->addWidget(m_resS, 8, 1);
    paramLLayout->addWidget(m_fibDoubTimeS, 9, 1);
    paramLLayout->addWidget(m_angS, 10, 1);
    paramLLayout->addWidget(m_angTimeS, 11, 1);
    paramLLayout->addWidget(m_DvegfS, 12, 1);
    paramLLayout->addWidget(m_VmaxVegfS, 13, 1);
    paramLLayout->addWidget(m_KmVegfS, 14, 1);
    paramLLayout->addWidget(m_hypVegfS, 15, 1);
    paramLLayout->addWidget(m_vegfThresS, 16, 1);

    QGridLayout *paramCLayout = new QGridLayout;
    paramCLayout->addWidget(m_treat, 1, 0);
    paramCLayout->addWidget(m_phaseRS, 2, 0);
    paramCLayout->addWidget(m_treatS, 1, 1);
    paramCLayout->addWidget(m_phaseRSS, 2, 1);
    paramCLayout->addWidget(m_alphaGroup, 3, 0, 1, 2);
    paramCLayout->addWidget(m_betaGroup, 4, 0, 1, 2);

    QGridLayout *paramRLayout = new QGridLayout;
    paramRLayout->addWidget(m_arrestTime, 1, 0);
    paramRLayout->addWidget(m_doseThres, 2, 0);
    paramRLayout->addWidget(m_fraction, 3, 0);
    paramRLayout->addWidget(m_oxy, 4, 0);
    paramRLayout->addWidget(m_D, 5, 0);
    paramRLayout->addWidget(m_Vmax, 6, 0);
    paramRLayout->addWidget(m_Km, 7, 0);
    paramRLayout->addWidget(m_pO2NormVes, 8, 0);
    paramRLayout->addWidget(m_pO2TumVes, 9, 0);
    paramRLayout->addWidget(m_hypThres, 10, 0);
    paramRLayout->addWidget(m_hypNecThres, 11, 0);

    paramRLayout->addWidget(m_arrestTimeS, 1, 1);
    paramRLayout->addWidget(m_doseThresS, 2, 1);
    paramRLayout->addWidget(m_fractionS, 3, 1);
    paramRLayout->addWidget(m_oxyS, 4, 1);
    paramRLayout->addWidget(m_DS, 5, 1);
    paramRLayout->addWidget(m_VmaxS, 6, 1);
    paramRLayout->addWidget(m_KmS, 7, 1);
    paramRLayout->addWidget(m_pO2NormVesS, 8, 1);
    paramRLayout->addWidget(m_pO2TumVesS, 9, 1);
    paramRLayout->addWidget(m_hypThresS, 10, 1);
    paramRLayout->addWidget(m_hypNecThresS, 11, 1);
    QHBoxLayout *paramLayout = new QHBoxLayout;
    paramLayout->addLayout(paramLLayout);
    paramLayout->addLayout(paramCLayout);
    paramLayout->addLayout(paramRLayout);
    m_param->setLayout(paramLayout);

    QHBoxLayout *hLayout = new QHBoxLayout;
    hLayout->addWidget(m_back);
    hLayout->addWidget(m_cancel);
    hLayout->addWidget(m_loadInData);
    hLayout->addWidget(m_start);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(m_param);
    mainLayout->addLayout(hLayout);
    mainLayout->addWidget(m_progress);

    setLayout(mainLayout);

    QObject::connect(m_back, SIGNAL(clicked()), this, SLOT(back()));
    QObject::connect(m_cancel, SIGNAL(clicked()), qApp, SLOT(quit()));
    QObject::connect(m_start, SIGNAL(clicked()), this, SLOT(start()));
    QObject::connect(m_loadInData, SIGNAL(clicked()), this, SLOT(selInDataFile()));

    loadInData("../InputFiles/refParMeanRT.dat");

    setWindowTitle("Radiotherapy Simulator");
    setWindowIcon(QIcon("../Figures/logo.png"));
    showMaximized();
}


void InAnWindow::back(){
    StartWindow *startWindow = new StartWindow;
    close();
    startWindow->show();
}


void InAnWindow::disable(int state){
}


int InAnWindow::loadInData(std::string nFInData){
    std::ifstream fInData(nFInData.c_str());

    double tumDens(0.0), sigmaTum(0.0), vascDens(0.0), sigmaVasc(0.0);

    //fInData >> tumDens >> sigmaTum >> vascDens >> sigmaVasc;

    m_tumDensS->setValue(tumDens);
    m_sigmaTumS->setValue(sigmaTum);
    m_vascDensS->setValue(vascDens);
    m_sigmaVascS->setValue(sigmaVasc);

    bool tumGrowth;
    int edgeOrder(0);
    double doubTime(0.0);

    fInData >> tumGrowth >> doubTime >> edgeOrder;

    m_doubTimeS->setValue(doubTime);
    m_edgeOrderS->setValue(edgeOrder);

    bool res;
    double fibDoubTime(0.0);

    fInData >> res;
    fInData >> fibDoubTime;
    m_fibDoubTimeS->setValue(fibDoubTime);

    bool ang;
    double angTime(0.0), vegfThres(0.0);

    fInData >> ang;
    fInData >> angTime >> vegfThres;

    m_angTimeS->setValue(angTime);
    m_vegfThresS->setValue(vegfThres);

    bool RT;
    int schedule(0);
    double fraction(0.0), interval(0.0), totalDose(0.0);
    double arrestTime(0.0), doseThres(0.0), hypNecThres(0.0);
    std::vector<double> alpha(8, 0.0);
    std::vector<double> beta(8, 0.0);

    fInData >> RT;
    for(int i(0); i < 8; i++){
        fInData >> alpha.at(i);
        fInData >> beta.at(i);
        beta.at(i) = alpha.at(i) / beta.at(i);
    }
    fInData >> doseThres  >> arrestTime >> hypNecThres;
    fInData >> fraction;

    for(int i(0); i < 8; i++){
        m_alphaS.at(i)->setValue(alpha.at(i));
        m_betaS.at(i)->setValue(beta.at(i));
    }
    m_doseThresS->setValue(doseThres);
    m_arrestTimeS->setValue(arrestTime);
    m_hypNecThresS->setValue(hypNecThres);
    m_fractionS->setValue(fraction);

    bool oxy;
    double constpO2(0.0), constpO2NormVes(0.0), constpO2TumVes(0.0);
    double D(0.0);
    double Km(0.0), pO2NormVes(0.0), pO2TumVes(0.0), Vmax(0.0);
    double Dvegf(0.0), VmaxVegf(0.0), KmVegf(0.0);
    double hypThres(0.0), hypVegf(0.0);

    fInData >> oxy;
    fInData >> D >> Vmax >> Km;
    fInData >> Dvegf >> VmaxVegf >> KmVegf;
    fInData >>  pO2NormVes >> pO2TumVes;
    fInData >> hypThres >> hypVegf;
    fInData.close();

    m_DS->setValue(D);
    m_VmaxS->setValue(Vmax);
    m_KmS->setValue(Km);
    m_DvegfS->setValue(Dvegf);
    m_VmaxVegfS->setValue(VmaxVegf);
    m_KmVegfS->setValue(KmVegf);
    m_pO2NormVesS->setValue(pO2NormVes);
    m_pO2TumVesS->setValue(pO2TumVes);
    m_hypThresS->setValue(hypThres);
    m_hypVegfS->setValue(hypVegf);

    return 0;
}


int InAnWindow::createInFiles(){
    QDir dir("../OutputFilesGUI");
    dir.removeRecursively();
    QDir().mkdir("../OutputFilesGUI");
    if(0){}
    /*if(m_histSpec->isChecked()){
        QFile::copy(QString::fromStdString("../HistSpec/tissueDim" +
                                           std::to_string(m_selHistSpec->currentIndex() + 1) +
                                           ".dat"), QString::fromStdString("../OutputFilesGUI/tissueDim.dat"));

        if (QFile::exists(QString::fromStdString("inTum.dat"))){
            QFile::remove(QString::fromStdString("inTum.dat"));
        }
        QFile::copy(QString::fromStdString("../HistSpec/inTum" +
                                           std::to_string(m_selHistSpec->currentIndex() + 1) +
                                           ".dat"), QString::fromStdString("inTum.dat"));

        if (QFile::exists(QString::fromStdString("inVes.dat"))){
            QFile::remove(QString::fromStdString("inVes.dat"));
        }
        QFile::copy(QString::fromStdString("../HistSpec/inVes" +
                                           std::to_string(m_selHistSpec->currentIndex() + 1) +
                                           ".dat"), QString::fromStdString("inVes.dat"));
    }*/

    else{

        int nrow(90), ncol(90), nlayer(1);
        double cellSize(20.0);
        std::ofstream fTissueDim("../OutputFilesGUI/tissueDim.dat");

        fTissueDim << nrow << std::endl;
        fTissueDim << ncol << std::endl;
        fTissueDim << nlayer   << std::endl;
        fTissueDim << cellSize << std::endl;

        fTissueDim.close();

        int ivd, ivd2;
        int mindim, mindim2, sqrtmin, tumToDist, vesToDist;
        int nrowNcol, nrowNcolNlayer;
        std::vector<int> div;
        std::vector<double> diff;

        nrowNcol = nrow * ncol;
        nrowNcolNlayer = nrowNcol * nlayer;
        tumToDist = m_tumDensS->value() * nrowNcolNlayer;

        if(m_vascDensS->value()){
            mindim = std::min(nrow, ncol);
            mindim2 = mindim * mindim;
            sqrtmin = sqrt(mindim) + 1;
            for(int l(1); l < sqrtmin; l++){
                if(!(nrow % l) && !(ncol % l)){
                    div.push_back(l);
                    diff.push_back(fabs(1.0 / (l * l) - m_vascDensS->value()));
                    div.push_back(mindim / l);
                    diff.push_back(fabs(double(l * l) / double(mindim2) - m_vascDensS->value()));
                }
            }

            ivd = div.at(std::min_element(diff.begin(), diff.end()) - diff.begin());
            ivd2 = ivd * ivd;
            vesToDist = std::min(nrowNcol / ivd2, nrowNcolNlayer - tumToDist);
        }

        else{
            vesToDist = 0;
        }

        int halfIvd, halfNrow, halfNcol, halfNlayer;
        int imHalfNrow2, lmHalfNlayer2;
        int lnrowNcol;
        int irr2, ishift;
        double R, RR;
        std::vector<simpCell> map(nrowNcolNlayer);

        halfIvd = 0.5 * ivd;
        halfNrow = 0.5 * nrow;
        halfNcol = 0.5 * ncol;
        halfNlayer = 0.5 * nlayer;
        R = sqrt(halfNrow * halfNrow + halfNcol * halfNcol + halfNlayer * halfNlayer);
        RR = 1.001 * sqrt(2.0 * halfIvd * halfIvd);

        int k(0);

        for(int l(0); l < nlayer; l++){
            lmHalfNlayer2 = (l - halfNlayer) * (l - halfNlayer);
            lnrowNcol = l * nrowNcol;
            for(int i(0); i < nrow; i++){
                imHalfNrow2 = (i - halfNrow) * (i - halfNrow);
                irr2 = (i % ivd - halfIvd) * (i % ivd - halfIvd);
                ishift = i / ivd * ncol / ivd;
                for(int j(0); j < ncol; j++){
                    map.at(k).r = sqrt(lmHalfNlayer2 + imHalfNrow2 +
                                       (j - halfNcol) * (j - halfNcol)) / R;
                    map.at(k).rr = lnrowNcol + ishift + j / ivd +
                            sqrt(irr2 + (j % ivd - halfIvd) * (j % ivd - halfIvd)) / RR;
                    map.at(k).k = k;
                    map.at(k).tum = 0;
                    map.at(k).ves = 0;
                    k++;
                }
            }
        }
        k = 0;

        std::sort(map.begin(), map.end(), compR);

        int m;
        double n;
        std::default_random_engine gen;
        std::normal_distribution<double> distTum(0, m_sigmaTumS->value());
        bool tooBigTumDens(m_tumDensS->value() > 0.9);
        bool tooSmallSigmaTum(m_sigmaTumS->value() < 0.15);

        if(!tooBigTumDens && !tooSmallSigmaTum){
            while(tumToDist > 0){
                n = distTum(gen);
                if(n >= 0.0 && n < 1.0){
                    m = n * nrowNcolNlayer;
                    if(!map.at(m).tum){
                        map.at(m).tum = 1;
                        tumToDist--;
                    }
                }
            }
        }

        else if(tooBigTumDens){
            for(int k(0); k < map.size(); k++){
                map.at(k).tum = 1;
            }
            std::reverse(map.begin(), map.end());
            int tumToRemove(nrowNcolNlayer - tumToDist);
            while(tumToRemove > 0){
                n = distTum(gen);
                if(n >= 0.0 && n < 1.0){
                    m = n * nrowNcolNlayer;
                    if(map.at(m).tum){
                        map.at(m).tum = 0;
                        tumToRemove--;
                    }
                }
            }
        }
        else{
            bool cond;
            while(tumToDist > 0){
                n = distTum(gen);
                if(n >= 0.0 && n < 1.0){
                    m = n * nrowNcolNlayer;
                    cond = true;
                    while(cond){
                        if(!map.at(m).tum){
                            map.at(m).tum = 1;
                            tumToDist--;
                            cond = false;
                        }
                        else{
                            m++;
                        }
                        if(m == nrowNcolNlayer){
                            m = 0;
                        }
                    }
                }
            }
        }

        std::sort(map.begin(), map.end(), compRR);

        int nivd2;
        bool availSpace, firstRound(true);
        std::normal_distribution<double> distVes(0, m_sigmaVascS->value());

        k = 0;
        while(vesToDist > 0){
            availSpace = false;
            for(int j(0); j < ivd2; j++){
                if(!map.at(k + j).ves && !map.at(k + j).tum){
                    availSpace = true;
                    break;
                }
            }

            if(availSpace){
                n = distVes(gen);
                while(n < 0.0 || n >= 1.0){
                    n = distVes(gen);
                }
                nivd2 = n * ivd2;
                m = k + nivd2;
                while(map.at(m).ves || map.at(m).tum){
                    nivd2++;
                    if(nivd2 >= ivd2){
                        nivd2 = 0;
                    }
                    m = k + nivd2;
                }
                map.at(m).ves = 1;
                vesToDist--;
                availSpace = false;
            }

            if(firstRound){
                k += ivd2;
            }
            else{
                k = rand() % (nrowNcol / ivd2) * ivd2;
            }
            if(k >= nrowNcol - ivd2){
                firstRound = false;
            }
        }

        std::sort(map.begin(), map.end(), compK);
        for(int k(0); k < nrowNcol; k++){
            if(map.at(k).ves == 1){
                for(int l(1); l < nlayer; l++){
                    map.at(l * nrowNcol + k).ves = 1;
                    map.at(l * nrowNcol + k).tum = 0;
                }
            }
        }

        std::ofstream fInTum("inTum.dat"), fInVes("inVes.dat");

        for(k = 0; k < nrowNcolNlayer; k++){
            fInTum << map.at(k).tum << std::endl;
            fInVes << map.at(k).ves << std::endl;
        }

        fInTum.close();
        fInVes.close();

    }

    std::ofstream fParam("../OutputFilesGUI/param.dat");
    std::vector<double> cycDistrib = {0.6, 0.25, 0.075, 0.075};
    for(int i(0); i < 4; i++){
        fParam << cycDistrib.at(i) << std::endl;
    }
//--------------------------------me quede aqui----------------------------------------
    /*fParam << m_paramTG->isChecked() << std::endl;
    if(m_paramTG->isChecked()){
        fParam << m_doubTime->value() << std::endl;
        fParam << m_edgeOrder->value() << std::endl;
        for(int i(0); i < 4; i++){
            fParam << m_time.at(i)->value() << std::endl;
        }
    }

    fParam << m_paramRes->isChecked() << std::endl;
    if(m_paramRes->isChecked()){
        fParam << m_fibDoubTime->value() << std::endl;
    }

    fParam << m_paramAng->isChecked() << std::endl;
    if(m_paramAng->isChecked()){
        fParam << m_angTime->value() << std::endl;
        fParam << m_Dvegf->value() *
                  m_oxySimTimeStep->value()  << std::endl;
        fParam << m_VmaxVegf->value() *
                  m_oxySimTimeStep-> value() << std::endl;
        fParam << m_KmVegf->value() << std::endl;
        fParam << m_hypVegf->value() << std::endl;
        fParam << m_vegfThres->value() << std::endl;
    }

    fParam << m_paramRT->isChecked() << std::endl;
    if(m_paramRT->isChecked()){
        for(int i(0); i < 8; i++){
            fParam << m_alpha.at(i)->value() << std::endl;
        }
        for(int i(0); i < 8; i++){
            fParam << m_beta.at(i)->value() << std::endl;
        }
        fParam << m_arrestTime->value() << std::endl;
        fParam << m_doseThres->value() << std::endl;
        fParam << m_fraction->value() << std::endl;
        fParam << m_totalDose->value() << std::endl;
        fParam << m_interval->value() << std::endl;
        if(m_MonFri->isChecked()){
            fParam << 0 << std::endl;
        }
        else{
            fParam << 1 << std::endl;
        }
    }

    fParam << m_paramOxy->isChecked() << std::endl;
    if(m_paramOxy->isChecked()){
        fParam << m_D->value() *
                  m_oxySimTimeStep->value() << std::endl;
        fParam << m_Vmax->value() *
                  m_oxySimTimeStep->value() << std::endl;
        fParam << m_Km->value() << std::endl;
        fParam << m_pO2NormVes->value() << std::endl;
        fParam << m_pO2TumVes->value() << std::endl;
        fParam << m_hypThres->value() << std::endl;
        fParam << m_hypNecThres->value() << std::endl;
    }
    else{
        fParam << m_constpO2->value() << std::endl;
        fParam << m_constpO2NormVes->value() << std::endl;
        fParam << m_constpO2TumVes->value() << std::endl;
    }

    fParam.close();

    std::ofstream fSimParam("../OutputFilesGUI/simParam.dat");
    m_simType = 1 * (m_nlayer->value() > 1) + 2 * m_coupSim;
    fSimParam << m_simType << std::endl;
    fSimParam << m_simTime->value() << std::endl;
    fSimParam << m_simTimeStep->value() << std::endl;
    fSimParam << m_oxySimTimeStep->value() << std::endl;

    fSimParam.close();*/
    return 0;
}


void InAnWindow::nextWindow(){
    m_progress->setRange(0, 0);
    new OutAnWindow("../OutputFilesGUI");

    close();
}


void InAnWindow::selInDataFile(){
    QString QnFInData = QFileDialog::getOpenFileName(this, "Select input data file", "../InputFiles");

    if(!QnFInData.isNull()){
        std::string nFInData = QnFInData.toUtf8().constData();
        loadInData(nFInData);
    }
}


void InAnWindow::start(){
    m_param->setEnabled(false);
    m_loadInData->setEnabled(false);
    m_start->setEnabled(false);
    //createInFiles();

    m_progress->reset();
    SimThread *thread = new SimThread(this);
    QObject::connect(thread, SIGNAL(resultReady(int)), this, SLOT(nextWindow(int)));
    QObject::connect(thread, SIGNAL(finished()), thread, SLOT(deleteLater()));
    QObject::connect(thread, SIGNAL(progressMax(int)), m_progress, SLOT(setMaximum(int)));
    QObject::connect(thread, SIGNAL(progress(int)), m_progress, SLOT(setValue(int)));
    thread->start();
}
