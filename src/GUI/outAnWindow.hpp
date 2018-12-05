#ifndef OUTANWINDOW_HPP
#define OUTANWINDOW_HPP

#include <QApplication>
#include <QDesktopWidget>
#include <QGroupBox>
#include <QImage>
#include <QLabel>
#include <QPixmap>
#include <QProgressBar>
#include <QPushButton>
#include <QStyle>
#include <QSpinBox>
#include <QtCharts>

#include <string>
#include <vector>

#include "inAnWindow.hpp"

class OutAnWindow : public QWidget{
    Q_OBJECT
public:
    OutAnWindow(std::string nFOutData);

private:
    int m_simNum;
    int m_simTime, m_simTimeStep;
    double m_endTreatTime, m_recTime;
    QComboBox *m_selChart;
    QGroupBox *m_selChartGroup;
    QChartView *m_chartView;
    QChart *m_cTumDens, *m_cTumVol, *m_cVascDens, *m_cKilledCells;
    QChart *m_cCycle, *m_cHypDens, *m_cPO2Stat, *m_cVegfStat;
    QValueAxis *m_xTumDens, *m_yTumDens;
    QValueAxis *m_xTumVol, *m_yTumVol;
    QValueAxis *m_xVascDens, *m_yVascDens;
    QValueAxis *m_xKilledCells, *m_yKilledCells;
    QValueAxis *m_xCycle, *m_yCycle;
    QValueAxis *m_xHypDens, *m_yHypDens;
    QValueAxis *m_xPO2Stat, *m_yPO2Stat;
    QValueAxis *m_xVegfStat, *m_yVegfStat;
    QLineSeries *m_endTreatDash, *m_recDash;
    QPushButton *m_sel, *m_change, *m_newSim;
    QPushButton *m_saveOutFiles;
    QPushButton *m_saveChart;
    QColor m_white, m_blueTum, m_red, m_redTum, m_black, m_blueTumDam;
    QColor m_yellowNec, m_greenApop, m_brownMit;
    QColor m_blueG1, m_green, m_orange, m_violet, m_redG0;

private slots:
    void change();
    void changeChart(const int numChart);
    void newSim();
    void saveChart();
    void saveOutFiles();
    void sel();
};

#endif
