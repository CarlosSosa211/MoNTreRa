/*#ifndef DEF_OUTWINDOW3D
#define DEF_OUTWINDOW3D

#include <Q3DScatter>
#include <QApplication>
#include <QDesktopWidget>
#include <QGroupBox>
#include <QLabel>
#include <QProgressBar>
#include <QPushButton>
#include <QStyle>
#include <QSpinBox>
#include <QtCharts>

#include <string>
#include <vector>

#include "inWindow.hpp"

class OutWindow3D : public QWidget{
    Q_OBJECT
public:
    OutWindow3D(std::string nFOutData);
    void drawChartDashLine(int iter);
    void drawMap(int numMap, int mapIter);

private:
    int m_mapNrow, m_mapNcol, m_mapNlayer, m_mapSclFac;
    int m_simTime, m_simTimeStep;
    std::vector<std::vector<int> > m_state, m_timer;
    std::vector<std::vector<double> > m_pO2;
    double m_maxpO2;
    QComboBox *m_selChart, *m_selMap;
    QGroupBox *m_selChartGroup, *m_selMapGroup, *m_mapGroup;
    QGroupBox *m_legendSt, *m_legendCyc, *m_legendOxy;
    QChartView *m_chartView;
    QChart *m_cTumDens, *m_cVascDens, *m_cKilledCells, *m_cCycle;
    QChart *m_cHypDens, *m_cPO2Stat, *m_cVegfStat;
    QValueAxis *m_xTumDens, *m_yTumDens;
    QValueAxis *m_xVascDens, *m_yVascDens;
    QValueAxis *m_xKilledCells, *m_yKilledCells;
    QValueAxis *m_xCycle, *m_yCycle;
    QValueAxis *m_xHypDens, *m_yHypDens;
    QValueAxis *m_xPO2Stat, *m_yPO2Stat;
    QValueAxis *m_xVegfStat, *m_yVegfStat;
    QLineSeries *m_sDash;
    QPushButton *m_sel, *m_change, *m_newSim;
    QPushButton *m_saveOutFiles;
    QPushButton *m_saveAllMaps, *m_saveChart, *m_saveMap;
    QLabel *m_time;
    QPushButton *m_play;
    bool m_stPlay;
    QSpinBox *m_timeS;
    QSlider *m_slider;
    QColor m_white, m_blueTum, m_red, m_redTum, m_black, m_blueTumDam;
    QColor m_yellowNec, m_greenApop, m_brownMit;
    QColor m_blueG1, m_green, m_orange, m_violet, m_redG0;
    QtDataVisualization::Q3DScatter *m_map;
    QWidget *m_mapCont;
    QtDataVisualization::QScatterDataProxy *m_proxyTum, *m_proxyTumDam, *m_proxyNormVes, *m_proxyTumVes;
    QtDataVisualization::QScatterDataProxy *m_proxyHypNec, *m_proxyMitCat, *m_proxyApop;
    QtDataVisualization::QScatterDataProxy *m_proxyG1, *m_proxyS, *m_proxyG2, *m_proxyM, *m_proxyG0;
    QtDataVisualization::QScatter3DSeries *m_sTum, *m_sTumDam, *m_sNormVes, *m_sTumVes;
    QtDataVisualization::QScatter3DSeries *m_sHypNec, *m_sMitCat, *m_sApop;
    QtDataVisualization::QScatter3DSeries *m_sG1, *m_sS, *m_sG2, *m_sM, *m_sG0;
    QtDataVisualization::QScatterDataArray *m_datTum, *m_datTumDam, *m_datNormVes, *m_datTumVes;
    QtDataVisualization::QScatterDataArray *m_datHypNec, *m_datMitCat, *m_datApop;
    QtDataVisualization::QScatterDataArray *m_datG1, *m_datS, *m_datG2, *m_datM, *m_datG0;
    QLabel *m_fib, *m_tum, *m_tumDam, *m_normVes, *m_tumVes;
    QLabel *m_hypNec, *m_mitCat, *m_apop;
    QLabel *m_fibSq, *m_tumSq, *m_tumDamSq, *m_normVesSq, *m_tumVesSq;
    QLabel *m_hypNecSq, *m_mitCatSq, *m_apopSq;
    QLabel *m_G1, *m_S, *m_G2, *m_M, *m_G0;
    QLabel *m_G1Sq, *m_SSq, *m_G2Sq, *m_MSq, *m_G0Sq;
    QLabel *m_oxyBar, *m_oxyValMax, *m_oxyValMin;

private slots:
    void change();
    void changeChart(int numChart);
    void changeIter(int iter);
    void changeNumMap(int numMap);
    void newSim();
    void play();
    void saveAllMaps();
    void saveChart();
    void saveMap();
    void saveOutFiles();
    void sel();

signals:
    void updateSlider(int mapIter);
};

#endif*/
