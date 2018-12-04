#ifndef DEF_STARTWINDOW
#define DEF_STARTWINDOW

#include <QApplication>
#include <QDesktopWidget>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>
#include <QRadioButton>
#include <QSlider>
#include <QStyle>

#include "inWindow.hpp"
#include "inAnWindow.hpp"

class StartWindow : public QWidget{
    Q_OBJECT
public:
    StartWindow();

private:
    QLabel *m_title;
    QGroupBox *m_sel;
    QRadioButton *m_fixedVal, *m_outFiles, *m_sensAn;
    QPushButton *m_ok;
    QProgressBar *m_progress;
    QLabel *m_progressOutL;

private slots:
    void advance();
};

#endif
