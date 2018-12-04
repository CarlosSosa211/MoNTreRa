#ifndef INANWINDOW_H
#define INANWINDOW_H

#include <string>

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QGroupBox>
#include <QProgressBar>
#include <QPushButton>
#include <QSpinBox>
#include <QWidget>

#include "startWindow.hpp"

class InAnWindow : public QWidget{
    Q_OBJECT
public:
    InAnWindow();
    int loadInData(std::string nFInData);

private:
    QGroupBox *m_param;
    QCheckBox *m_tumGrowth, *m_res, *m_ang, *m_treat, *m_phaseRS, *m_oxy;
    QCheckBox *m_tumDens, *m_sigmaTum, *m_vascDens, *m_sigmaVasc;
    QCheckBox *m_doubTime;
    QCheckBox *m_edgeOrder;
    QGroupBox *m_timeGroup, *m_distGroup;
    QCheckBox *m_fibDoubTime;
    QCheckBox *m_angTime, *m_Dvegf, *m_vegfThres;
    QCheckBox *m_VmaxVegf, *m_KmVegf, *m_hypVegf;
    QGroupBox *m_alphaGroup, *m_betaGroup;
    std::vector<QCheckBox *> m_alpha, m_beta;
    QCheckBox *m_arrestTime, *m_doseThres;
    QCheckBox *m_fraction;
    QCheckBox *m_D, *m_Vmax, *m_Km, *m_pO2NormVes, *m_pO2TumVes;
    QCheckBox *m_hypThres, *m_hypNecThres;
    QDoubleSpinBox *m_tumDensS, *m_sigmaTumS, *m_vascDensS, *m_sigmaVascS;
    QDoubleSpinBox *m_doubTimeS;
    QSpinBox *m_edgeOrderS;
    QDoubleSpinBox *m_fibDoubTimeS;
    QDoubleSpinBox *m_angTimeS, *m_DvegfS, *m_vegfThresS;
    QDoubleSpinBox *m_VmaxVegfS, *m_KmVegfS, *m_hypVegfS;
    std::vector<QDoubleSpinBox *> m_alphaS, m_betaS;
    QDoubleSpinBox *m_arrestTimeS, *m_doseThresS;
    QDoubleSpinBox *m_fractionS;
    QDoubleSpinBox *m_DS, *m_VmaxS, *m_KmS, *m_pO2NormVesS, *m_pO2TumVesS;
    QDoubleSpinBox *m_hypThresS, *m_hypNecThresS;
    QPushButton *m_loadInData;
    QPushButton *m_back, *m_cancel, *m_start;
    QProgressBar *m_progress;

private slots:
    void back();
    void start();
    void disable(int state);
    void selInDataFile();
};

#endif
