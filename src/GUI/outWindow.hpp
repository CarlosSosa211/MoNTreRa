#ifndef DEF_OUTWINDOW
#define DEF_OUTWINDOW

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

#include "inWindow.hpp"

class OutWindow : public QWidget{
    Q_OBJECT
public:
    OutWindow(std::string nFOutData);
    void drawChartDashLine(const int iter);
    void drawMap(const int numMap, const int mapIter);

private:
    int m_mapNrow, m_mapNcol, m_mapSclFac;
    int m_simTime, m_simTimeStep;
    double m_endTreatTime, m_rec, m_recTime;
    std::vector<std::vector<int> > m_state, m_timer;
    std::vector<std::vector<double> > m_pO2, m_vegf;
    double m_maxpO2, m_maxvegf;
    QComboBox *m_selChart, *m_selMap;
    QGroupBox *m_selChartGroup, *m_selMapGroup, *m_mapGroup;
    QGroupBox *m_legendSt, *m_legendCyc, *m_legendPO2, *m_legendVegf;
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
    QLineSeries *m_sDash, *m_endTreatDash, *m_recDash;
    QPushButton *m_sel, *m_change, *m_newSim;
    QPushButton *m_saveOutFiles;
    QPushButton *m_saveAllMaps, *m_saveChart, *m_saveMap;
    QLabel *m_map, *m_time;
    QImage *m_imState, *m_imTimer, *m_imPO2, *m_imVegf;
    QPixmap *m_pixState, *m_pixTimer, *m_pixPO2, *m_pixVegf;
    QPushButton *m_play;
    bool m_stPlay;
    QSpinBox *m_timeS;
    QSlider *m_slider;
    QColor m_white, m_blueTum, m_red, m_redTum, m_black, m_blueTumDam;
    QColor m_yellowNec, m_greenApop, m_brownMit;
    QColor m_blueG1, m_green, m_orange, m_violet, m_redG0;
    QLabel *m_fib, *m_tum, *m_tumDam, *m_normVes, *m_tumVes, *m_hypNec;
    QLabel *m_mitCat, *m_apop;
    QLabel *m_fibSq, *m_tumSq, *m_tumDamSq, *m_normVesSq, *m_tumVesSq;
    QLabel *m_hypNecSq, *m_mitCatSq, *m_apopSq;
    QLabel *m_G1, *m_S, *m_G2, *m_M, *m_G0;
    QLabel *m_G1Sq, *m_SSq, *m_G2Sq, *m_MSq, *m_G0Sq;
    QLabel *m_pO2Bar, *m_pO2ValMax, *m_pO2ValMin;
    QLabel *m_vegfBar, *m_vegfValMax, *m_vegfValMin;

private slots:
    void change();
    void changeChart(const int numChart);
    void changeIter(const int iter);
    void changeNumMap(const int numMap);
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

#endif

