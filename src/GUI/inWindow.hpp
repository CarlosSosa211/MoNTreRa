#ifndef DEF_INWINDOW
#define DEF_INWINDOW

#include <string>
#include <vector>

#include <QApplication>
#include <QComboBox>
#include <QDesktopWidget>
#include <QDoubleSpinBox>
#include <QGroupBox>
#include <QLabel>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QRadioButton>
#include <QSpinBox>
#include <QStyle>

#include "../constOxyTissue.hpp"
#include "../coupler.hpp"
#include "outWindow.hpp"
#include "outWindow3D.hpp"
#include "outWindowOxy.hpp"
#include "../oxyTissue.hpp"
#include "../rootSimulator.hpp"
#include "../simpcell.hpp"
#include "../simulator.hpp"
#include "startWindow.hpp"
#include "../treatment.hpp"
#include "../tissue.hpp"
#include "simThread.hpp"

class InWindow : public QWidget{
    Q_OBJECT
public:
    InWindow();
    int createInFiles();
    int loadInData(std::string nFInData);

public slots:
    int simulate();

private:
    QGroupBox *m_param, *m_paramArch, *m_paramTG, *m_paramRes;
    QGroupBox *m_paramAng, *m_paramRT, *m_paramRS, *m_paramTreat;
    QGroupBox *m_paramOxy, *m_paramSim;
    QRadioButton *m_artif, *m_histSpec;
    QComboBox *m_selHistSpec;
    QGroupBox *m_artifGroup;
    QSpinBox *m_nrow, *m_ncol, *m_nlayer;
    QDoubleSpinBox *m_cellSize;
    QDoubleSpinBox *m_tumDens, *m_sigmaTum, *m_vascDens, *m_sigmaVasc;
    QDoubleSpinBox *m_doubTime;
    QSpinBox *m_edgeOrder;
    QGroupBox *m_timeGroup, *m_distGroup;
    std::vector<QDoubleSpinBox *> m_time, m_dist;
    QDoubleSpinBox *m_fibDoubTime;
    QDoubleSpinBox *m_angTime, *m_Dvegf, *m_vegfThres;
    QDoubleSpinBox *m_VmaxVegf, *m_KmVegf, *m_hypVegf;
    QGroupBox *m_alphaGroup, *m_betaGroup;
    std::vector<QDoubleSpinBox *> m_alpha, m_beta;
    QDoubleSpinBox *m_arrestTime, *m_doseThres;
    QDoubleSpinBox *m_fraction, *m_totalDose, *m_interval;
    QGroupBox *m_schedule;
    QRadioButton *m_MonFri, *m_everyDay;
    QDoubleSpinBox *m_DO2, *m_VmaxO2, *m_KmO2, *m_pO2NormVes, *m_pO2TumVes;
    QDoubleSpinBox *m_hypThres, *m_hypNecThres;
    QDoubleSpinBox *m_constpO2, *m_constpO2NormVes, *m_constpO2TumVes;
    QLabel *m_simTimeL, *m_oxySimTimeL, *m_simTimeStepL, *m_oxySimTimeStepL;
    QSpinBox *m_simTime, *m_simTimeStep, *m_oxySimTimeStep;
    QPushButton *m_loadInData;
    QPushButton *m_back, *m_cancel, *m_simulate;
    QProgressBar *m_progress;
    QLabel *m_progressCreaL, *m_progressSimL, *m_progressOutL;
    int m_simType;
    bool m_coupSim;

private slots:
    void back();
    void changeSimTime(bool value);
    void changeSimTime(double value);
    void coupSim(bool value);
    void disable(bool oxy);
    void nextWindow(int simType);
    void selInDataFile();

};

#endif
