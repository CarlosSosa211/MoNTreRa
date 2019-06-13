#include <fstream>
#include <stdlib.h>
#include <time.h>

#include <QGridLayout>
#include <QVBoxLayout>

#include "startWindow.hpp"

/*------------------------------------------------------------------------------
 * Constructor of the class StartWindow.
------------------------------------------------------------------------------*/

StartWindow::StartWindow() : QWidget(){
    m_title        = new QLabel;
    m_sel          = new QGroupBox("Select mode", this);
    m_fixedVal     = new QRadioButton("Simulation with fixed values", m_sel);
    m_outFiles     = new QRadioButton("Load output files", m_sel);
    m_sensAn       = new QRadioButton("Sensitivity analysis", m_sel);
    m_ok           = new QPushButton("Ok", this);
    m_progress     = new QProgressBar(this);
    m_progressOutL = new QLabel("Preparing outputs...", this);

    m_title->setPixmap(QPixmap("../Figures/lLogo.png"));
    m_fixedVal->setChecked(true);
    m_progressOutL->hide();

    QVBoxLayout *boxLayout  = new QVBoxLayout;
    QGridLayout *mainLayout = new QGridLayout;

    boxLayout->addWidget(m_fixedVal);
    boxLayout->addWidget(m_outFiles);
    boxLayout->addWidget(m_sensAn);
    m_sel->setLayout(boxLayout);

    mainLayout->addWidget(m_title, 0, 0, 1, 2, Qt::AlignHCenter);
    mainLayout->addWidget(m_sel, 2, 0, 1, 2);
    mainLayout->addWidget(m_ok, 3, 0, 1, 2);
    mainLayout->addWidget(m_progress, 4, 0, 1, 1);
    mainLayout->addWidget(m_progressOutL, 4, 1, 1, 1);
    setLayout(mainLayout);

    QObject::connect(m_ok, SIGNAL(clicked()), this, SLOT(advance()));
    srand(time(NULL));

    setWindowTitle("Radiotherapy Simulator");
    setWindowIcon(QIcon("../Figures/logo.png"));
    resize(500, 300);
    setGeometry(QStyle::alignedRect(Qt::LeftToRight, Qt::AlignCenter, size(),
                                    qApp->desktop()->availableGeometry()));
    show();
}


/*------------------------------------------------------------------------------
 * This slot goes to the corresponding window (InWindow, OutWindowOxy,
 * OutWindow or OutWindow3D).
------------------------------------------------------------------------------*/

void StartWindow::advance(){
    if(m_fixedVal->isChecked()){
        new InWindow;
        close();
    }

    else if(m_outFiles->isChecked()){
        QString QnFOutData = QFileDialog::getExistingDirectory(
                    this, "Select output data folder",
                    "../", QFileDialog::ShowDirsOnly |
                    QFileDialog::DontResolveSymlinks);

        if(!QnFOutData.isNull()){
            m_progress->setRange(0, 0);
            m_progressOutL->show();
            m_sel->setEnabled(false);
            m_ok->setEnabled(false);

            std::string nFOutData = QnFOutData.toUtf8().constData();
            std::ifstream fSimParam((nFOutData + "/simParam.dat").c_str());
            int simType;

            fSimParam >> simType;
            fSimParam.close();

            switch(simType){
            case 0:{
                new OutWindowOxy(nFOutData);
                break;
            }

            case 2:{
                new OutWindow(nFOutData);
                break;
            }

            case 3:{
                //new OutWindow3D(nFOutData);
                break;
            }
            }
            close();
        }
    }

    else if(m_sensAn->isChecked()){
        new InAnWindow;
        close();
    }
}

