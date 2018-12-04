#include <fstream>
#include <iostream>

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

    //m_tumGrowthS = new QDoubleSpinBox(m_param);
    m_doubTimeS  = new QDoubleSpinBox(m_param);
    m_edgeOrderS = new QSpinBox(m_param);

    //m_resS = new QDoubleSpinBox(m_param);
    m_fibDoubTimeS = new QDoubleSpinBox( m_param);

    //m_angS       = new QDoubleSpinBox(m_param);
    m_angTimeS   = new QDoubleSpinBox(m_param);
    m_DvegfS     = new QDoubleSpinBox(m_param);
    m_VmaxVegfS  = new QDoubleSpinBox(m_param);
    m_KmVegfS    = new QDoubleSpinBox(m_param);
    m_hypVegfS   = new QDoubleSpinBox(m_param);
    m_vegfThresS = new QDoubleSpinBox(m_param);

    //m_treatS      = new QDoubleSpinBox(m_param);
    //m_phaseRS    = new QDoubleSpinBox(m_param);
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

    //m_oxyS         = new QDoubleSpinBox(m_param);
    m_DS           = new QDoubleSpinBox(m_param);
    m_VmaxS        = new QDoubleSpinBox(m_param);
    m_KmS          = new QDoubleSpinBox(m_param);
    m_pO2NormVesS  = new QDoubleSpinBox(m_param);
    m_pO2TumVesS   = new QDoubleSpinBox(m_param);
    m_hypThresS    = new QDoubleSpinBox(m_param);
    m_hypNecThresS = new QDoubleSpinBox(m_param);

    m_back     = new QPushButton("Back", this);
    m_cancel   = new QPushButton("Cancel", this);
    m_loadInData = new QPushButton("Load input data", this);
    m_start    = new QPushButton("Start analysis", this);
    m_progress = new QProgressBar(this);

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
    //paramLLayout->addWidget(m_tumGrowthS, 5, 0);
    paramLLayout->addWidget(m_doubTimeS, 6, 1);
    paramLLayout->addWidget(m_edgeOrderS, 7, 1);
    //paramLLayout->addWidget(m_resS, 8, 0);
    paramLLayout->addWidget(m_fibDoubTimeS, 9, 1);
    //paramLLayout->addWidget(m_angS, 10, 0);
    paramLLayout->addWidget(m_angTimeS, 11, 1);
    paramLLayout->addWidget(m_DvegfS, 12, 1);
    paramLLayout->addWidget(m_VmaxVegfS, 13, 1);
    paramLLayout->addWidget(m_KmVegfS, 14, 1);
    paramLLayout->addWidget(m_hypVegfS, 15, 1);
    paramLLayout->addWidget(m_vegfThresS, 16, 1);


    QGridLayout *paramCLayout = new QGridLayout;
    paramCLayout->addWidget(m_treat, 1, 0);
    paramCLayout->addWidget(m_phaseRS, 2, 0);
    paramCLayout->addWidget(m_alphaGroup, 3, 0);
    paramCLayout->addWidget(m_betaGroup, 4, 0);

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
    //paramRLayout->addWidget(m_oxyS, 4, 1);
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
    }
    for(int i(0); i < 8; i++){
        fInData >> beta.at(i);
    }
    fInData >> doseThres  >> arrestTime >> hypNecThres;
    fInData >> fraction;

    for(int i(0); i < 8; i++){
        m_alphaS.at(i)->setValue(alpha.at(i));
        m_betaS.at(i)->setValue(beta.at(i));
    }
    m_arrestTimeS->setValue(arrestTime);
    m_doseThresS->setValue(doseThres);
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
    m_pO2NormVesS->setValue(pO2NormVes);
    m_pO2TumVesS->setValue(pO2TumVes);
    m_hypThresS->setValue(hypThres);
    m_hypNecThresS->setValue(hypNecThres);

    return 0;
}


void InAnWindow::selInDataFile(){
    QString QnFInData = QFileDialog::getOpenFileName(this, "Select input data file", "../InputFiles");

    if(!QnFInData.isNull()){
        std::string nFInData = QnFInData.toUtf8().constData();
        loadInData(nFInData);
    }
}


void InAnWindow::start(){
}
