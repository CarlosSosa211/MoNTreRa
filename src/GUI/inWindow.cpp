#include <algorithm>
#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <stdlib.h>
#include <string>
#include <time.h>

#include <QFormLayout>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>

#include "inWindow.hpp"

InWindow::InWindow() : QWidget(){
    m_param      = new QGroupBox("Select the processes to be simulated and enter the value of the parameters", this);
    m_paramArch  = new QGroupBox("Tissue architecture");
    m_paramTG    = new QGroupBox("Tumor growth", m_param);
    m_paramRes   = new QGroupBox("Fibroblasts proliferation", m_param);
    m_paramAng   = new QGroupBox("Angiogenesis", m_param);
    m_paramRT    = new QGroupBox("Response to radiotherapy", m_param);
    m_paramRS    = new QGroupBox("Radiosensitivity", m_paramRT);
    m_paramTreat = new QGroupBox("Treatment", m_paramRT);
    m_paramOxy   = new QGroupBox("Oxygenation", m_param);
    m_paramSim   = new QGroupBox("Enter the simulation parameters", this);

    m_histSpec    = new QRadioButton("Histopathological specimen", m_paramArch);
    m_artif       = new QRadioButton("Artificial tissue ", m_paramArch);
    m_selHistSpec = new QComboBox(m_histSpec);
    m_artifGroup  = new QGroupBox(m_paramArch);
    m_nrow        = new QSpinBox(m_artifGroup);
    m_ncol        = new QSpinBox(m_artifGroup);
    m_nlayer      = new QSpinBox(m_artifGroup);
    m_tumDens     = new QDoubleSpinBox(m_artifGroup);
    m_sigmaTum    = new QDoubleSpinBox(m_artifGroup);
    m_vascDens    = new QDoubleSpinBox(m_artifGroup);
    m_sigmaVasc   = new QDoubleSpinBox(m_paramArch);
    m_distGroup   = new QGroupBox("Initial distribution of cells", m_paramArch);
    for(int i(0); i < 4; i++){
        m_dist.push_back(new QDoubleSpinBox(m_distGroup));
    }

    m_paramTG->setCheckable(true);
    m_doubTime  = new QDoubleSpinBox(m_paramTG);
    m_edgeOrder  = new QSpinBox(m_paramTG);
    m_timeGroup = new QGroupBox("Fraction of duration of each phase", m_paramTG);
    for(int i(0); i < 4; i++){
        m_time.push_back(new QDoubleSpinBox(m_timeGroup));
    }

    m_paramRes->setCheckable(true);
    m_fibDoubTime = new QDoubleSpinBox(m_paramRes);

    m_paramAng->setCheckable(true);
    m_angTime   = new QDoubleSpinBox(m_paramAng);
    m_Dvegf     = new QDoubleSpinBox(m_paramAng);
    m_VmaxVegf  = new QDoubleSpinBox(m_paramAng);
    m_KmVegf    = new QDoubleSpinBox(m_paramAng);
    m_hypVegf   = new QDoubleSpinBox(m_paramAng);
    m_vegfThres = new QDoubleSpinBox(m_paramAng);


    m_paramRT->setCheckable(true);
    m_alphaGroup = new QGroupBox("Alpha (Gy⁻¹)", m_paramRS);
    m_betaGroup  = new QGroupBox("Beta (Gy⁻²)", m_paramRS);
    for(int i(0); i < 8; i++){
        m_alpha.push_back(new QDoubleSpinBox(m_alphaGroup));
        m_beta.push_back(new QDoubleSpinBox(m_betaGroup));
    }
    m_arrestTime = new QDoubleSpinBox(m_paramRS);

    m_fraction  = new QDoubleSpinBox(m_paramTreat);
    m_totalDose = new QDoubleSpinBox(m_paramTreat);
    m_interval  = new QDoubleSpinBox(m_paramTreat);
    m_schedule  = new QGroupBox("Schedule", m_paramTreat);
    m_MonFri    = new QRadioButton("Mon - Fri", m_schedule);
    m_everyDay  = new QRadioButton("Every day", m_schedule);

    m_paramOxy->setCheckable(true);
    m_D           = new QDoubleSpinBox(m_paramOxy);
    m_Vmax        = new QDoubleSpinBox(m_paramOxy);
    m_Km          = new QDoubleSpinBox(m_paramOxy);
    m_pO2NormVes  = new QDoubleSpinBox(m_paramOxy);
    m_pO2TumVes   = new QDoubleSpinBox(m_paramOxy);
    m_hypThres    = new QDoubleSpinBox(m_paramOxy);
    m_hypNecThres = new QDoubleSpinBox(m_paramOxy);

    m_constpO2        = new QDoubleSpinBox(m_param);
    m_constpO2NormVes = new QDoubleSpinBox(m_param);
    m_constpO2TumVes  = new QDoubleSpinBox(m_param);

    m_simTimeL        = new QLabel("Simulation time (h)", m_paramSim);
    m_oxySimTimeL     = new QLabel("Simulation time (ms)", m_paramSim);
    m_simTimeStepL    = new QLabel("Simulation timestep (h)", m_paramSim);
    m_oxySimTimeStepL = new QLabel("Oxygenation simulation timestep (ms)", m_paramSim);

    m_simTime        = new QSpinBox(m_paramSim);
    m_simTimeStep    = new QSpinBox(m_paramSim);
    m_oxySimTimeStep = new QSpinBox(m_paramSim);

    m_back       = new QPushButton("Back", this);
    m_cancel     = new QPushButton("Cancel", this);
    m_loadInData = new QPushButton("Load input data", this);
    m_simulate   = new QPushButton("Simulate", this);

    m_progress      = new QProgressBar(this);
    m_progressCreaL = new QLabel("Generating initial tissue...");
    m_progressSimL  = new QLabel("Simulating...", this);
    m_progressOutL  = new QLabel("Preparing outputs...", this);

    m_histSpec->setChecked(true);

    QStringList histSpecList = {"Tissue P14 01685 6-1-cd31 101x62",
                                "Tissue P14 01685 7-5-cd31 97x57",
                                "Tissue P14 01685 8-5-cd31 111x75",
                                "Tissue P14 08021 7-2-cd31 125x73",
                                "Tissue P14 08021 8-2-cd31 120x71",
                                "Tissue P14 08021 9-2-cd31 91x55",
                                "Tissue P14 12769 12-2-cd31 97x56",
                                "Tissue P14 12769 13-4-cd31 108x65",
                                "Tissue P14 12769 14-2-cd31 109x64",
                                "Tissue P14 20984 5-1-cd31 95x58",
                                "Tissue P14 20984 6-1-cd31 119x70",
                                "Tissue P14 20984 7-1-cd31 107x66",
                                "Tissue P14 26519 5-3-cd31 93x55",
                                "Tissue P14 26519 6-3-cd31 55x37",
                                "Tissue P14 26519 7-3-cd31 95x59",
                                "Tissue P14 26671 4-1-cd31 94x58",
                                "Tissue P14 26671 5-3-cd31 104x65",
                                "Tissue P14 26671 6-1-cd31 98x63",
                                "Tissue P14 29860 4-1-cd31 101x56",
                                "Tissue P14 29860 5-1-cd31 113x64",
                                "Tissue P14 29860 6-1-cd31 125x71"};
    m_selHistSpec->addItems(histSpecList);
    m_artif->setChecked(false);
    m_artifGroup->setEnabled(false);
    m_nrow->setMaximum(499);
    m_ncol->setMaximum(499);
    m_nlayer->setMaximum(499);
    m_tumDens->setMaximum(1.0);
    m_tumDens->setSingleStep(0.01);
    m_sigmaTum->setMaximum(2.0);
    m_sigmaTum->setSingleStep(0.01);
    m_vascDens->setMaximum(1.0);
    m_vascDens->setSingleStep(0.01);
    m_sigmaVasc->setMaximum(2.0);
    m_sigmaVasc->setSingleStep(0.01);
    for(int i(0); i < 4; i++){
        m_dist.at(i)->setMaximum(1.0);
        m_dist.at(i)->setDecimals(3);
        m_dist.at(i)->setSingleStep(0.01);
    }

    m_doubTime->setMaximum(9999.9);
    m_doubTime->setSingleStep(1.0);
    m_doubTime->setDecimals(1);
    m_edgeOrder->setMaximum(10);
    for(int i(0); i < 4; i++){
        m_time.at(i)->setMaximum(1.0);
        m_time.at(i)->setDecimals(3);
        m_time.at(i)->setSingleStep(0.01);
    }

    m_fibDoubTime->setMaximum(9999.9);
    m_fibDoubTime->setSingleStep(1.0);
    m_fibDoubTime->setDecimals(1);

    m_angTime->setMaximum(9999.9);
    m_angTime->setSingleStep(1.0);
    m_angTime->setDecimals(1);
    m_Dvegf->setMaximum(9999.9);
    m_Dvegf->setSingleStep(0.1);
    m_Dvegf->setDecimals(1);
    m_VmaxVegf->setDecimals(3);
    m_KmVegf->setDecimals(3);
    m_hypVegf->setDecimals(1);
    m_vegfThres->setDecimals(1);

    for(int i(0); i < 8; i++){
        m_alpha.at(i)->setSingleStep(0.001);
        m_alpha.at(i)->setDecimals(3);
        m_beta.at(i)->setSingleStep(0.001);
        m_beta.at(i)->setDecimals(3);
    }

    m_fraction->setDecimals(1);
    m_totalDose->setDecimals(1);
    m_interval->setDecimals(1);

    m_D->setMaximum(9999.9);
    m_D->setSingleStep(0.1);
    m_D->setDecimals(1);
    m_Vmax->setDecimals(3);
    m_Km->setDecimals(3);
    m_hypNecThres->setDecimals(3);

    m_constpO2->setEnabled(false);
    m_constpO2NormVes->setEnabled(false);
    m_constpO2TumVes->setEnabled(false);

    m_simTime->setMaximum(20000);

    QFormLayout *formLayoutArch     = new QFormLayout;
    QFormLayout *formLayoutTG       = new QFormLayout;
    QFormLayout *formLayoutRes      = new QFormLayout;
    QFormLayout *formLayoutAng      = new QFormLayout;
    QFormLayout *formLayoutRS       = new QFormLayout;
    QFormLayout *formLayoutOxy      = new QFormLayout;
    QFormLayout *formLayoutTreat    = new QFormLayout;
    QFormLayout *formLayoutConstOxy = new QFormLayout;

    QFormLayout *artifLayout = new QFormLayout;
    artifLayout->addRow("Number of rows", m_nrow);
    artifLayout->addRow("Number of columns", m_ncol);
    artifLayout->addRow("Number of layers", m_nlayer);
    artifLayout->addRow("Initial tumor density", m_tumDens);
    artifLayout->addRow("Standard deviation\nof tumor cells", m_sigmaTum);
    artifLayout->addRow("Initial vascular density", m_vascDens);
    artifLayout->addRow("Standard deviation\nof vessels", m_sigmaVasc);
    m_artifGroup->setLayout(artifLayout);

    QFormLayout *distLayout = new QFormLayout;
    distLayout->addRow("Phase G1", m_dist.at(0));
    distLayout->addRow("Phase S", m_dist.at(1));
    distLayout->addRow("Phase G2", m_dist.at(2));
    distLayout->addRow("Phase M", m_dist.at(3));
    m_distGroup->setLayout(distLayout);

    formLayoutArch->addRow(m_histSpec);
    formLayoutArch->addRow(m_selHistSpec);
    formLayoutArch->addRow(m_artif);
    formLayoutArch->addRow(m_artifGroup);
    formLayoutArch->addRow(m_distGroup);
    m_paramArch->setLayout(formLayoutArch);

    formLayoutTG->addRow("Doubling time of tumor cells (h)", m_doubTime);
    formLayoutTG->addRow("Edge size", m_edgeOrder);
    QFormLayout *timeLayout = new QFormLayout;
    timeLayout->addRow("Phase G1", m_time.at(0));
    timeLayout->addRow("Phase S", m_time.at(1));
    timeLayout->addRow("Phase G2", m_time.at(2));
    timeLayout->addRow("Phase M", m_time.at(3));
    m_timeGroup->setLayout(timeLayout);
    formLayoutTG->addRow(m_timeGroup);
    m_paramTG->setLayout(formLayoutTG);

    formLayoutRes->addRow("Doubling time of fibroblasts (h)", m_fibDoubTime);
    m_paramRes->setLayout(formLayoutRes);

    formLayoutAng->addRow("Doubling time of vessels (h)", m_angTime);
    formLayoutAng->addRow("Diffusion coefficient (μm²/ms)", m_Dvegf);
    formLayoutAng->addRow("Vmax (mol/μm³ms)", m_VmaxVegf);
    formLayoutAng->addRow("K (mol/μm³)", m_KmVegf);
    formLayoutAng->addRow("Hypoxic cells VEGF (mol/μm³)", m_hypVegf);
    formLayoutAng->addRow("VEGF threshold (mol/μm³)", m_vegfThres);
    m_paramAng->setLayout(formLayoutAng);

    QFormLayout *alphaLayout = new QFormLayout;
    alphaLayout->addRow("Fibroblasts", m_alpha.at(0));
    alphaLayout->addRow("Tumor cells in phase G1", m_alpha.at(1));
    alphaLayout->addRow("Tumor cells in phase S", m_alpha.at(2));
    alphaLayout->addRow("Tumor cells in phase G2", m_alpha.at(3));
    alphaLayout->addRow("Tumor cells in phase M", m_alpha.at(4));
    alphaLayout->addRow("Tumor cells in phase G0", m_alpha.at(5));
    alphaLayout->addRow("Normal vessels", m_alpha.at(6));
    alphaLayout->addRow("Tumor vessels", m_alpha.at(7));
    m_alphaGroup->setLayout(alphaLayout);
    formLayoutRS->addRow(m_alphaGroup);

    QFormLayout *betaLayout = new QFormLayout;
    betaLayout->addRow("Fibroblasts", m_beta.at(0));
    betaLayout->addRow("Tumor cells in phase G1", m_beta.at(1));
    betaLayout->addRow("Tumor cells in phase S", m_beta.at(2));
    betaLayout->addRow("Tumor cells in phase G2", m_beta.at(3));
    betaLayout->addRow("Tumor cells in phase M", m_beta.at(4));
    betaLayout->addRow("Tumor cells in phase G0", m_beta.at(5));
    betaLayout->addRow("Normal vessels", m_beta.at(6));
    betaLayout->addRow("Tumor vessels", m_beta.at(7));
    m_betaGroup->setLayout(betaLayout);
    formLayoutRS->addRow(m_betaGroup);
    formLayoutRS->addRow("Time of cycle arrest (h)", m_arrestTime);

    m_paramRS->setLayout(formLayoutRS);

    formLayoutTreat->addRow("Fraction (Gy)", m_fraction);
    formLayoutTreat->addRow("Total dose (Gy)", m_totalDose);
    formLayoutTreat->addRow("Interval (h)", m_interval);
    QVBoxLayout *schLayout = new QVBoxLayout;
    schLayout->addWidget(m_MonFri);
    schLayout->addWidget(m_everyDay);
    m_schedule->setLayout(schLayout);
    formLayoutTreat->addRow(m_schedule);
    m_paramTreat->setLayout(formLayoutTreat);

    QHBoxLayout *layoutRT = new QHBoxLayout;
    layoutRT->addWidget(m_paramRS);
    layoutRT->addWidget(m_paramTreat);
    m_paramRT->setLayout(layoutRT);

    formLayoutOxy->addRow("Diffusion coefficient (μm²/ms)", m_D);
    formLayoutOxy->addRow("Vmax (mmHg/ms)", m_Vmax);
    formLayoutOxy->addRow("K (mmHg)", m_Km);
    formLayoutOxy->addRow("Normal vessels pO2 (mmHg)", m_pO2NormVes);
    formLayoutOxy->addRow("Tumor vessels pO2 (mmHg)", m_pO2TumVes);
    formLayoutOxy->addRow("Hypoxia threshold (mmHg)", m_hypThres);
    formLayoutOxy->addRow("Hypoxic necrosis threshold (mmHg)", m_hypNecThres);
    m_paramOxy->setLayout(formLayoutOxy);

    formLayoutConstOxy->addRow("Constant pO2 (mmHg)", m_constpO2);
    formLayoutConstOxy->addRow("Normal vessels constant pO2 (mmHg)", m_constpO2NormVes);
    formLayoutConstOxy->addRow("Tumor vessels constant pO2 (mmHg)", m_constpO2TumVes);

    QVBoxLayout *vOxyLayout = new QVBoxLayout;
    vOxyLayout->addWidget(m_paramOxy);
    vOxyLayout->addLayout(formLayoutConstOxy);

    QGridLayout *layoutSim = new QGridLayout;
    layoutSim->addWidget(m_simTimeL, 0, 0);
    layoutSim->addWidget(m_oxySimTimeL, 0, 0);
    layoutSim->addWidget(m_simTime, 0, 1);
    layoutSim->addWidget(m_simTimeStepL, 0, 2);
    layoutSim->addWidget(m_simTimeStep, 0, 3);
    layoutSim->addWidget(m_oxySimTimeStepL, 0, 4);
    layoutSim->addWidget(m_oxySimTimeStep, 0, 5);
    m_paramSim->setLayout(layoutSim);

    QVBoxLayout *vboxLayout = new QVBoxLayout;
    vboxLayout->addWidget(m_paramTG);
    vboxLayout->addWidget(m_paramRes);
    vboxLayout->addWidget(m_paramAng);

    QHBoxLayout *hboxLayout = new QHBoxLayout;
    hboxLayout->addWidget(m_paramArch);
    hboxLayout->addLayout(vboxLayout);
    hboxLayout->addWidget(m_paramRT);
    hboxLayout->addLayout(vOxyLayout);
    m_param->setLayout(hboxLayout);

    QHBoxLayout *hLayout = new QHBoxLayout;
    hLayout->addWidget(m_back);
    hLayout->addWidget(m_cancel);
    hLayout->addWidget(m_loadInData);
    hLayout->addWidget(m_simulate);

    QGridLayout *progressLayout = new QGridLayout;
    progressLayout->addWidget(m_progress, 0, 0);
    progressLayout->addWidget(m_progressCreaL, 0, 1);
    progressLayout->addWidget(m_progressSimL, 0, 1);
    progressLayout->addWidget(m_progressOutL, 0, 1);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(m_param);
    mainLayout->addWidget(m_paramSim);
    mainLayout->addLayout(hLayout);
    mainLayout->addLayout(progressLayout);

    setLayout(mainLayout);

    m_progressCreaL->hide();
    m_progressSimL->hide();
    m_progressOutL->hide();

    QObject::connect(m_histSpec, SIGNAL(toggled(bool)), m_selHistSpec, SLOT(setEnabled(bool)));
    QObject::connect(m_artif, SIGNAL(toggled(bool)), m_artifGroup, SLOT(setEnabled(bool)));
    QObject::connect(m_fraction, SIGNAL(valueChanged(double)), this, SLOT(changeSimTime(double)));
    QObject::connect(m_totalDose, SIGNAL(valueChanged(double)), this, SLOT(changeSimTime(double)));
    QObject::connect(m_interval, SIGNAL(valueChanged(double)), this, SLOT(changeSimTime(double)));
    QObject::connect(m_MonFri, SIGNAL(clicked(bool)), this, SLOT(changeSimTime(bool)));
    QObject::connect(m_everyDay, SIGNAL(clicked(bool)), this, SLOT(changeSimTime(bool)));
    QObject::connect(m_paramTG, SIGNAL(toggled(bool)), this, SLOT(coupSim(bool)));
    QObject::connect(m_paramRes, SIGNAL(toggled(bool)), this, SLOT(coupSim(bool)));
    QObject::connect(m_paramAng, SIGNAL(toggled(bool)), this, SLOT(coupSim(bool)));
    QObject::connect(m_paramRT, SIGNAL(toggled(bool)), this, SLOT(coupSim(bool)));
    QObject::connect(m_paramOxy, SIGNAL(toggled(bool)), this, SLOT(disable(bool)));
    QObject::connect(m_loadInData, SIGNAL(clicked()), this, SLOT(selInDataFile()));
    QObject::connect(m_back, SIGNAL(clicked()), this, SLOT(back()));
    QObject::connect(m_cancel, SIGNAL(clicked()), qApp, SLOT(quit()));
    QObject::connect(m_simulate, SIGNAL(clicked()), this, SLOT(simulate()));

    loadInData("../InputFiles/inTG.dat");

    setWindowTitle("Radiotherapy Simulator");
    setWindowIcon(QIcon("../Figures/logo.png"));
    showMaximized();
}


void InWindow::back(){
    new StartWindow;
    close();
}


void InWindow::changeSimTime(bool value){
    changeSimTime(0.0);
}


void InWindow::changeSimTime(double value){
    if(m_MonFri->isChecked()){
        m_simTime->setValue((int(m_totalDose->value() / m_fraction->value()) / 5 * 7.0 +
                             int(m_totalDose->value() / m_fraction->value()) % 5) *
                            m_interval->value());
    }
    if(m_everyDay->isChecked()){
        m_simTime->setValue((m_totalDose->value() / m_fraction->value()) *
                            m_interval->value());
    }
}


void InWindow::coupSim(bool value){
    m_coupSim = m_paramTG->isChecked() || m_paramRes->isChecked() ||
            m_paramRT->isChecked();
    if(m_coupSim){
        m_simTimeL->show();
        m_oxySimTimeL->hide();
    }

    else{
        m_simTimeL->hide();
        m_oxySimTimeL->show();
    }

    m_simTimeStep->setEnabled(m_coupSim);
}


int InWindow::createInFiles(){
    QDir dir("../OutputFilesGUI");
    dir.removeRecursively();
    QDir().mkdir("../OutputFilesGUI");
    if(m_histSpec->isChecked()){
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
    }

    else{
        std::ofstream fTissueDim("../OutputFilesGUI/tissueDim.dat");

        fTissueDim << m_nrow->value() << std::endl;
        fTissueDim << m_ncol->value() << std::endl;
        fTissueDim << m_nlayer->value() << std::endl;

        fTissueDim.close();

        int ivd, ivd2;
        int mindim, mindim2, sqrtmin, tumToDist, vesToDist;
        int nrowNcol, nrowNcolNlayer;
        std::vector<int> div;
        std::vector<double> diff;

        nrowNcol = m_nrow->value() * m_ncol->value();
        nrowNcolNlayer = nrowNcol * m_nlayer->value();
        tumToDist = m_tumDens->value() * nrowNcolNlayer;

        if(m_vascDens->value()){
            mindim = std::min(m_nrow->value(), m_ncol->value());
            mindim2 = mindim * mindim;
            sqrtmin = sqrt(mindim) + 1;
            for(int l(1); l < sqrtmin; l++){
                if(!(m_nrow->value() % l) && !(m_ncol->value() % l)){
                    div.push_back(l);
                    diff.push_back(fabs(1.0 / (l * l) - m_vascDens->value()));
                    div.push_back(mindim / l);
                    diff.push_back(fabs(double(l * l) / double(mindim2) - m_vascDens->value()));
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
        halfNrow = 0.5 * m_nrow->value();
        halfNcol = 0.5 * m_ncol->value();
        halfNlayer = 0.5 * m_nlayer->value();
        R = sqrt(halfNrow * halfNrow + halfNcol * halfNcol + halfNlayer * halfNlayer);
        RR = 1.001 * sqrt(2.0 * halfIvd * halfIvd);

        int k(0);

        for(int l(0); l < m_nlayer->value(); l++){
            lmHalfNlayer2 = (l - halfNlayer) * (l - halfNlayer);
            lnrowNcol = l * nrowNcol;
            for(int i(0); i < m_nrow->value(); i++){
                imHalfNrow2 = (i - halfNrow) * (i - halfNrow);
                irr2 = (i % ivd - halfIvd) * (i % ivd - halfIvd);
                ishift = i / ivd * m_ncol->value() / ivd;
                for(int j(0); j < m_ncol->value(); j++){
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
        std::normal_distribution<double> distTum(0, m_sigmaTum->value());
        bool tooBigTumDens(m_tumDens->value() > 0.9);
        bool tooSmallSigmaTum(m_sigmaTum->value() < 0.15);

        if(!tooBigTumDens && !tooSmallSigmaTum){
            m_progress->setMaximum(vesToDist + tumToDist);
            while(tumToDist > 0){
                n = distTum(gen);
                if(n >= 0.0 && n < 1.0){
                    m = n * nrowNcolNlayer;
                    if(!map.at(m).tum){
                        map.at(m).tum = 1;
                        tumToDist--;
                        m_progress->setValue(m_progress->maximum() - tumToDist - vesToDist);
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
            m_progress->setMaximum(vesToDist + tumToRemove);
            while(tumToRemove > 0){
                n = distTum(gen);
                if(n >= 0.0 && n < 1.0){
                    m = n * nrowNcolNlayer;
                    if(map.at(m).tum){
                        map.at(m).tum = 0;
                        tumToRemove--;
                        m_progress->setValue(m_progress->maximum() - tumToRemove - vesToDist);
                    }
                }
            }
        }
        else{
            m_progress->setMaximum(vesToDist + tumToDist);
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
                            m_progress->setValue(m_progress->maximum() - tumToDist - vesToDist);
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
        std::normal_distribution<double> distVes(0, m_sigmaVasc->value());

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
                m_progress->setValue(m_progress->maximum() - vesToDist);
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
                for(int l(1); l < m_nlayer->value(); l++){
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

        m_progress->setValue(m_progress->maximum());
    }

    std::ofstream fParam("../OutputFilesGUI/param.dat");

    for(int i(0); i < 4; i++){
        fParam << m_dist.at(i)->value() << std::endl;
    }

    fParam << m_paramTG->isChecked() << std::endl;
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

    double cellSize(20.0); //um
    double nu(m_oxySimTimeStep->value() / (cellSize * cellSize));

    fParam << m_paramAng->isChecked() << std::endl;
    if(m_paramAng->isChecked()){
        fParam << m_angTime->value() << std::endl;
        fParam << m_Dvegf->value() * nu << std::endl;
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
        fParam << m_D->value() * nu << std::endl;
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

    fSimParam.close();
    return 0;
}


void InWindow::disable(bool oxy){
    m_constpO2NormVes->setEnabled(!oxy);
    m_constpO2TumVes->setEnabled(!oxy);
    m_constpO2->setEnabled(!oxy);
}


int InWindow::loadInData(std::string nFInData){
    m_simTimeL->hide();
    m_oxySimTimeL->show();
    m_simTimeStep->setEnabled(false);

    std::ifstream fInData(nFInData.c_str());

    bool histSpec;
    int nrow(0), ncol(0), nlayer(0);
    double tumDens(0.0), sigmaTum(0.0), vascDens(0.0), sigmaVasc(0.0);

    fInData >> histSpec;
    m_histSpec->setChecked(histSpec);
    m_artif->setChecked(!histSpec);
    if(histSpec){
        m_selHistSpec->setCurrentIndex(histSpec - 1);
    }
    else{
        fInData >> nrow >> ncol >> nlayer;
        fInData >> tumDens >> sigmaTum >> vascDens >> sigmaVasc;
    }

    m_nrow->setValue(nrow);
    m_ncol->setValue(ncol);
    m_nlayer->setValue(nlayer);
    m_tumDens->setValue(tumDens);
    m_sigmaTum->setValue(sigmaTum);
    m_vascDens->setValue(vascDens);
    m_sigmaVasc->setValue(sigmaVasc);

    std::vector<double> cycDistrib(4, 0.0);

    for(int i(0); i < 4; i++){
        fInData >> cycDistrib.at(i);
    }
    for(int i(0); i < 4; i++){
        m_dist.at(i)->setValue(cycDistrib.at(i));
    }

    bool tumGrowth;
    int edgeOrder(0);
    double doubTime(0.0);
    std::vector<double> cycDur(4, 0.0);

    fInData >> tumGrowth;
    m_paramTG->setChecked(tumGrowth);
    if(tumGrowth){
        fInData >> doubTime >> edgeOrder;
        for(int i(0); i < 4; i++){
            fInData >> cycDur.at(i);
        }
        m_simTimeL->show();
        m_oxySimTimeL->hide();
        m_simTimeStep->setEnabled(true);
    }

    m_doubTime->setValue(doubTime);
    m_edgeOrder->setValue(edgeOrder);
    for(int i(0); i < 4; i++){
        m_time.at(i)->setValue(cycDur.at(i));
    }

    bool res;
    double fibDoubTime(0.0);

    fInData >> res;
    m_paramRes->setChecked(res);
    if(res){
        fInData >> fibDoubTime;
        m_simTimeL->show();
        m_oxySimTimeL->hide();
        m_simTimeStep->setEnabled(true);
    }
    m_fibDoubTime->setValue(fibDoubTime);

    bool ang;
    double angTime(0.0), Dvegf(0.0), VmaxVegf(0.0), KmVegf(0.0);
    double hypVegf(0.0), vegfThres(0.0);

    fInData >> ang;
    m_paramAng->setChecked(ang);
    if(ang){
        fInData >> angTime >> Dvegf >> VmaxVegf >> KmVegf;
        fInData >> hypVegf >> vegfThres;
        m_simTimeL->show();
        m_oxySimTimeL->hide();
        m_simTimeStep->setEnabled(true);
    }

    m_angTime->setValue(angTime);
    m_Dvegf->setValue(Dvegf);
    m_VmaxVegf->setValue(VmaxVegf);
    m_KmVegf->setValue(KmVegf);
    m_hypVegf->setValue(hypVegf);
    m_vegfThres->setValue(vegfThres);

    bool RT;
    int schedule(0);
    double fraction(0.0), interval(0.0), totalDose(0.0);
    double arrestTime(0.0);
    std::vector<double> alpha(8, 0.0);
    std::vector<double> beta(8, 0.0);

    fInData >> RT;
    m_paramRT->setChecked(RT);
    if(RT){
        for(int i(0); i < 8; i++){
            fInData >> alpha.at(i);
        }
        for(int i(0); i < 8; i++){
            fInData >> beta.at(i);
        }
        fInData >> arrestTime;
        fInData >> fraction >> totalDose >> interval >> schedule;
        m_simTimeL->show();
        m_oxySimTimeL->hide();
        m_simTimeStep->setEnabled(true);
    }

    for(int i(0); i < 8; i++){
        m_alpha.at(i)->setValue(alpha.at(i));
        m_beta.at(i)->setValue(beta.at(i));
    }
    m_arrestTime->setValue(arrestTime);
    m_fraction->setValue(fraction);
    m_totalDose->setValue(totalDose);
    m_interval->setValue(interval);
    if(!schedule){
        m_MonFri->setChecked(true);
    }
    else{
        m_everyDay->setChecked(true);
    }

    bool oxy;
    double constpO2(0.0), constpO2NormVes(0.0), constpO2TumVes(0.0);
    double D(0.0), hypThres(0.0), hypNecThres(0.0);
    double Km(0.0), pO2NormVes(0.0), pO2TumVes(0.0), Vmax(0.0);

    fInData >> oxy;
    m_paramOxy->setChecked(oxy);
    m_constpO2->setEnabled(!oxy);
    m_constpO2NormVes->setEnabled(!oxy);
    m_constpO2TumVes->setEnabled(!oxy);
    if(oxy){
        fInData >> D >> Vmax >> Km >> pO2NormVes;
        fInData >> pO2TumVes >> hypThres >> hypNecThres;
    }
    else{
        fInData >> constpO2 >> constpO2NormVes >> constpO2TumVes;
    }

    m_D->setValue(D);
    m_Vmax->setValue(Vmax);
    m_Km->setValue(Km);
    m_pO2NormVes->setValue(pO2NormVes);
    m_pO2TumVes->setValue(pO2TumVes);
    m_hypThres->setValue(hypThres);
    m_hypNecThres->setValue(hypNecThres);

    m_constpO2->setValue(constpO2);
    m_constpO2NormVes->setValue(constpO2NormVes);
    m_constpO2TumVes->setValue(constpO2TumVes);

    int oxySimTimeStep, simTime, simTimeStep;

    fInData >> simTime >> simTimeStep >> oxySimTimeStep;
    m_simTime->setValue(simTime);
    m_simTimeStep->setValue(simTimeStep);
    m_oxySimTimeStep->setValue(oxySimTimeStep);

    fInData.close();

    m_coupSim = m_paramTG->isChecked() || m_paramRes->isChecked() ||
            m_paramRT->isChecked();

    return 0;
}


void InWindow::nextWindow(int simType){
    m_progress->setRange(0, 0);
    m_progressSimL->hide();
    m_progressOutL->show();

    switch(simType){
    case 0:{
        new OutWindowOxy("../OutputFilesGUI");
        break;
    }

    case 2:{
        new OutWindow("../OutputFilesGUI");
        break;
    }

    case 3:{
        //new OutWindow3D("../OutputFilesGUI");
        break;
    }
    }
    close();
}


void InWindow::selInDataFile(){
    QString QnFInData = QFileDialog::getOpenFileName(this, "Select input data file", "../InputFiles");

    if(!QnFInData.isNull()){
        std::string nFInData = QnFInData.toUtf8().constData();
        loadInData(nFInData);
    }
}


int InWindow::simulate(){
    m_paramArch->setEnabled(false);
    m_paramTG->setEnabled(false);
    m_paramRes->setEnabled(false);
    m_paramAng->setEnabled(false);
    m_paramRT->setEnabled(false);
    m_paramOxy->setEnabled(false);
    m_constpO2->setEnabled(false);
    m_constpO2NormVes->setEnabled(false);
    m_constpO2TumVes->setEnabled(false);
    m_paramSim->setEnabled(false);
    m_loadInData->setEnabled(false);
    m_simulate->setEnabled(false);

    m_progressCreaL->show();
    createInFiles();
    m_progressCreaL->hide();
    m_progressSimL->show();
    m_progress->reset();
    SimThread *thread = new SimThread(this);
    QObject::connect(thread, SIGNAL(resultReady(int)), this, SLOT(nextWindow(int)));
    QObject::connect(thread, SIGNAL(finished()), thread, SLOT(deleteLater()));
    QObject::connect(thread, SIGNAL(progressMax(int)), m_progress, SLOT(setMaximum(int)));
    QObject::connect(thread, SIGNAL(progress(int)), m_progress, SLOT(setValue(int)));
    thread->start();
    return 0;
}
