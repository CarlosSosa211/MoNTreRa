#ifndef DEF_OUTWINDOWOXY
#define DEF_OUTWINDOWOXY

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

class OutWindowOxy : public QWidget{
    Q_OBJECT
public:
    OutWindowOxy(std::string nFOutData);
    void drawChartDashLine(int iter);
    void drawMap(int numMap, int mapIter);

private:
    int m_mapNrow, m_mapNcol, m_mapSclFac;
    int m_simTime, m_simOxyTimeStep;
    double m_stableTime;
    std::vector<std::vector<double> > m_pO2, m_vegf;
    double m_maxpO2, m_maxvegf;
    QComboBox *m_selChart, *m_selMap;
    QGroupBox *m_selChartGroup, *m_selMapGroup, *m_mapGroup;
    QGroupBox *m_legendPO2, *m_legendVegf;
    QChartView *m_chartView;
    QChart *m_cHypDens, *m_cPO2Stat, *m_cVegfStat;
    QValueAxis *m_xHypDens, *m_yHypDens;
    QValueAxis *m_xPO2Stat, *m_yPO2Stat;
    QValueAxis *m_xVegfStat, *m_yVegfStat;
    QLineSeries *m_sDash, *m_stableDash;
    QPushButton *m_sel, *m_change, *m_newSim;
    QPushButton *m_saveOutFiles;
    QPushButton *m_saveAllMaps, *m_saveChart, *m_saveMap;
    QLabel *m_map, *m_time;
    QImage *m_imPO2, *m_imVegf;
    QPixmap *m_pixPO2, *m_pixVegf;
    QPushButton *m_play;
    bool m_stPlay;
    QSpinBox *m_timeS;
    QSlider *m_slider;
    QColor m_green;
    QLabel *m_pO2Bar, *m_pO2ValMax, *m_pO2ValMin;
    QLabel *m_vegfBar, *m_vegfValMax, *m_vegfValMin;

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

#endif

