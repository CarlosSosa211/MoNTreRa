#include <fstream>
#include <iostream>

#include <QGridLayout>
#include <QHBoxLayout>
#include <QMessageBox>
#include <QVBoxLayout>

#include "outWindowOxy.hpp"

OutWindowOxy::OutWindowOxy(std::string nFOutData) : QWidget(){
    m_selChartGroup = new QGroupBox("Select a chart", this);
    m_selMapGroup = new QGroupBox("Select a map", this);
    m_mapGroup = new QGroupBox(this);

    m_selChart = new QComboBox(m_selChartGroup);
    m_selMap = new QComboBox(m_selMapGroup);

    m_chartView = new QChartView;

    double a, b, c;

    std::ifstream fHypDens((nFOutData + "/hypDens.res").c_str());

    if(!fHypDens.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening hypDens.res");
    }

    else{
        m_cHypDens = new QChart;
        m_xHypDens = new QValueAxis;
        m_yHypDens = new QValueAxis;
        m_cHypDens->addAxis(m_xHypDens, Qt::AlignBottom);
        m_cHypDens->addAxis(m_yHypDens, Qt::AlignLeft);
        m_cHypDens->setMaximumWidth(1000);

        QLineSeries *sHypDens = new QLineSeries(m_cHypDens);

        fHypDens >> a >> b;
        while(!fHypDens.eof()){
            sHypDens->append(a, b);
            fHypDens >> a >> b;
        }
        fHypDens.close();

        m_cHypDens->addSeries(sHypDens);
        m_cHypDens->legend()->hide();
        m_cHypDens->setTitle("Evolution of the hypoxic density");

        m_xHypDens->setTitleText("Time (ms)");
        m_xHypDens->setLabelFormat("%i");
        sHypDens->attachAxis(m_xHypDens);
        m_xHypDens->applyNiceNumbers();

        m_yHypDens->setTitleText("Hypoxic density");
        m_yHypDens->setLabelFormat("%i");
        sHypDens->attachAxis(m_yHypDens);
        m_yHypDens->setMin(0.0);
        m_yHypDens->applyNiceNumbers();

        m_sDash = new QLineSeries(m_cHypDens);
        m_sDash->append(0.0, 0.0);
        m_sDash->append(0.0, 100.0);
        m_sDash->setPen(QPen(Qt::DashLine));
        m_cHypDens->addSeries(m_sDash);
        m_sDash->attachAxis(m_xHypDens);
        m_sDash->attachAxis(m_yHypDens);

        m_selChart->addItem("Hypoxic density");
    }

    std::ifstream fPO2Stat((nFOutData + "/pO2Stat.res").c_str());

    if(!fPO2Stat.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening pO2Stat.res");
    }

    else{
        m_cPO2Stat = new QChart;
        m_xPO2Stat = new QValueAxis;
        m_yPO2Stat = new QValueAxis;
        m_cPO2Stat->addAxis(m_xPO2Stat, Qt::AlignBottom);
        m_cPO2Stat->addAxis(m_yPO2Stat, Qt::AlignLeft);
        m_cPO2Stat->setMaximumWidth(1000);

        QLineSeries *sPO2Med = new QLineSeries(m_cPO2Stat);
        QLineSeries *sPO2Mean = new QLineSeries(m_cPO2Stat);

        fPO2Stat >> a >> b >> c;
        while(!fPO2Stat.eof()){
            sPO2Med->append(a, b);
            sPO2Mean->append(a, c);
            fPO2Stat >> a >> b >> c;
        }
        fPO2Stat.close();

        m_cPO2Stat->addSeries(sPO2Med);
        m_cPO2Stat->addSeries(sPO2Mean);
        m_cPO2Stat->setTitle("Evolution of the PO2 statistics");
        sPO2Med->setName("Median");
        sPO2Mean->setName("Mean");

        m_xPO2Stat->setTitleText("Time (ms)");
        m_xPO2Stat->setLabelFormat("%i");
        sPO2Mean->attachAxis(m_xPO2Stat);
        sPO2Med->attachAxis(m_xPO2Stat);
        m_xPO2Stat->applyNiceNumbers();

        m_yPO2Stat->setTitleText("Value");
        m_yPO2Stat->setLabelFormat("%4.2f");
        sPO2Mean->attachAxis(m_yPO2Stat);
        sPO2Med->attachAxis(m_yPO2Stat);
        m_yPO2Stat->setMin(0.0);
        m_yPO2Stat->applyNiceNumbers();

        m_selChart->addItem("PO2 statistics");
    }

    std::ifstream fVegfStat((nFOutData + "/vegfStat.res").c_str());

    if(!fVegfStat.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening vegfStat.res");
    }

    else{
        m_cVegfStat = new QChart;
        m_xVegfStat = new QValueAxis;
        m_yVegfStat = new QValueAxis;
        m_cVegfStat->addAxis(m_xVegfStat, Qt::AlignBottom);
        m_cVegfStat->addAxis(m_yVegfStat, Qt::AlignLeft);
        m_cVegfStat->setMaximumWidth(1000);

        QLineSeries *sVegfMed = new QLineSeries(m_cVegfStat);
        QLineSeries *sVegfMean = new QLineSeries(m_cVegfStat);

        fVegfStat >> a >> b >> c;
        while(!fVegfStat.eof()){
            sVegfMed->append(a, b);
            sVegfMean->append(a, c);
            fVegfStat >> a >> b >> c;
        }
        fVegfStat.close();

        m_cVegfStat->addSeries(sVegfMed);
        m_cVegfStat->addSeries(sVegfMean);
        m_cVegfStat->setTitle("Evolution of the VEGF statistics");
        sVegfMed->setName("Median");
        sVegfMean->setName("Mean");

        m_xVegfStat->setTitleText("Time (h)");
        m_xVegfStat->setLabelFormat("%i");
        sVegfMean->attachAxis(m_xVegfStat);
        sVegfMed->attachAxis(m_xVegfStat);
        m_xVegfStat->applyNiceNumbers();

        m_yVegfStat->setTitleText("Value");
        m_yVegfStat->setLabelFormat("%4.2f");
        sVegfMean->attachAxis(m_yVegfStat);
        sVegfMed->attachAxis(m_yVegfStat);
        m_yVegfStat->setMin(0.0);
        m_yVegfStat->applyNiceNumbers();

        m_selChart->addItem("VEGF statistics");
    }

    m_chartView->setChart(m_cHypDens);

    m_map = new QLabel(m_mapGroup);

    std::ifstream fTissueDim((nFOutData + "/tissueDim.dat").c_str());

    if(!fTissueDim.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening tissueDim.dat");
    }

    else{
        fTissueDim >> m_mapNrow >> m_mapNcol;
        fTissueDim.close();
    }

    int nrowNcol(m_mapNrow * m_mapNcol);

    m_mapSclFac = 16;

    std::ifstream fSimParam((nFOutData + "/simParam.dat").c_str());

    if(!fSimParam.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening simParam.dat");
    }

    else{
        fSimParam >> m_simTime >> m_simTime >> m_simOxyTimeStep >> m_simOxyTimeStep;
        fSimParam.close();
    }

    std::ifstream fPO2((nFOutData + "/po2.res").c_str());

    if(!fPO2.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening po2.res");
    }

    else{
        int i(0);
        double tempd;

        m_maxpO2 = 0.0;

        fPO2 >> tempd;
        while(!fPO2.eof()){
            m_pO2.push_back(std::vector<double>());
            for(int j(0); j < nrowNcol; j++){
                m_pO2[i].push_back(tempd);
                if(tempd > m_maxpO2){
                    m_maxpO2 = tempd;
                }
                fPO2 >> tempd;
            }
            qApp->processEvents();
            i++;
        }
        fPO2.close();

        m_pixPO2 = new QPixmap;

        m_legendPO2 = new QGroupBox(m_mapGroup);

        m_imPO2 = new QImage(m_mapSclFac * m_mapNcol, m_mapSclFac * m_mapNrow,
                             QImage::Format_RGB32);
        m_pO2Bar = new QLabel(m_legendPO2);
        m_pO2ValMax = new QLabel(QString::number(m_maxpO2) + " mmHg");
        m_pO2ValMin = new QLabel("0 mmHg");

        int nrowPO2(50), ncolPO2(720);
        QImage pO2Im = QImage(ncolPO2, nrowPO2, QImage::Format_RGB32);

        for(int i(0); i < nrowPO2; i++){
            for(int j(0); j < ncolPO2; j++){
                pO2Im.setPixelColor(j, i, QColor::fromHsv((j/3) % 240, 200, 255));
            }
        }
        m_pO2Bar->setPixmap(QPixmap::fromImage(pO2Im));

        QGridLayout *legOxyLayout = new QGridLayout();
        legOxyLayout->addWidget(m_pO2Bar, 0, 1, 1, 2);
        legOxyLayout->addWidget(m_pO2ValMax, 0, 0);
        legOxyLayout->addWidget(m_pO2ValMin, 0, 3);

        m_legendPO2->setMinimumWidth(100);
        m_legendPO2->setLayout(legOxyLayout);

        m_selMap->addItem("pO2");
    }


    std::ifstream fVegf((nFOutData + "/vegf.res").c_str());

    if(!fVegf.is_open()){
        QMessageBox::critical(this, "Error", "Problem while opening vegf.res");
    }

    else{
        int i(0);
        double tempd;

        m_maxvegf = 0.0;

        fVegf >> tempd;
        while(!fVegf.eof()){
            m_vegf.push_back(std::vector<double>());
            for(int j(0); j < nrowNcol; j++){
                m_vegf[i].push_back(tempd);
                if(tempd > m_maxvegf){
                    m_maxvegf = tempd;
                }
                fVegf >> tempd;
            }
            qApp->processEvents();
            i++;
        }
        fVegf.close();

        m_pixVegf = new QPixmap;

        m_legendVegf = new QGroupBox(m_mapGroup);

        m_imVegf = new QImage(m_mapSclFac * m_mapNcol, m_mapSclFac * m_mapNrow,
                              QImage::Format_RGB32);
        m_vegfBar = new QLabel(m_legendVegf);
        m_vegfValMax = new QLabel(QString::number(m_maxvegf) + "mol/μm³");
        m_vegfValMin = new QLabel("0 mol/μm³");

        int nrowVegf(50), ncolVegf(720);
        QImage vegfIm = QImage(ncolVegf, nrowVegf, QImage::Format_RGB32);

        for(int i(0); i < nrowVegf; i++){
            for(int j(0); j < ncolVegf; j++){
                vegfIm.setPixelColor(j, i, QColor::fromHsv((j/3) % 240, 200, 255));
            }
        }
        m_vegfBar->setPixmap(QPixmap::fromImage(vegfIm));

        QGridLayout *legVegfLayout = new QGridLayout();
        legVegfLayout->addWidget(m_vegfBar, 0, 1, 1, 2);
        legVegfLayout->addWidget(m_vegfValMax, 0, 0);
        legVegfLayout->addWidget(m_vegfValMin, 0, 3);

        m_legendVegf->setMinimumWidth(100);
        m_legendVegf->setLayout(legVegfLayout);

        m_selMap->addItem("VEGF");
    }

    m_time   = new QLabel("Time (ms)", m_mapGroup);
    m_timeS  = new QSpinBox(m_mapGroup);
    m_play   = new QPushButton(m_mapGroup);
    m_slider = new QSlider(Qt::Horizontal, this);

    m_timeS->setMinimum(0);
    m_timeS->setMaximum((m_pO2.size() - 1) * m_simOxyTimeStep);
    m_slider->setMinimum(0);
    m_slider->setMaximum((m_pO2.size() - 1) * m_simOxyTimeStep);
    m_play->setIcon(QIcon("../Figures/play.png"));
    m_play->setFixedSize(40, 40);
    m_play->setIconSize(QSize(30, 30));
    m_stPlay = false;

    m_sel          = new QPushButton("Select a new mode", this);
    m_change       = new QPushButton("Change input and parameters", this);
    m_saveOutFiles = new QPushButton("Save output files", this);
    m_saveChart    = new QPushButton("Save chart", this);
    m_saveMap      = new QPushButton("Save map of current iteration", this);
    m_saveAllMaps  = new QPushButton("Save maps of every iteration", this);
    m_newSim       = new QPushButton("New simulation", this);

    QVBoxLayout *chartLayout = new QVBoxLayout;
    chartLayout->addWidget(m_selChart);
    m_selChartGroup->setLayout(chartLayout);

    QVBoxLayout *selMapLayout = new QVBoxLayout;
    selMapLayout->addWidget(m_selMap);
    m_selMapGroup->setLayout(selMapLayout);


    QHBoxLayout *timeLayout = new QHBoxLayout;
    timeLayout->addWidget(m_time);
    timeLayout->addWidget(m_timeS);

    QGridLayout *mapLayout = new QGridLayout;
    mapLayout->addWidget(m_map, 0, 0, 1, 2);
    mapLayout->addWidget(m_legendPO2, 1, 0, 1, 2);
    mapLayout->addWidget(m_legendVegf, 1, 0, 1, 2);
    mapLayout->addWidget(m_slider, 2, 0);
    mapLayout->addLayout(timeLayout, 3, 0, 1, 2);
    mapLayout->addWidget(m_play, 2, 1, 1, 1, Qt::AlignCenter);
    m_mapGroup->setLayout(mapLayout);

    QHBoxLayout *hLayout = new QHBoxLayout;
    hLayout->addWidget(m_sel);
    hLayout->addWidget(m_change);
    hLayout->addWidget(m_saveOutFiles);
    hLayout->addWidget(m_saveChart);
    hLayout->addWidget(m_saveMap);
    hLayout->addWidget(m_saveAllMaps);
    hLayout->addWidget(m_newSim);

    QGridLayout *gridLayout = new QGridLayout;
    gridLayout->addWidget(m_selChartGroup, 0, 0);
    gridLayout->addWidget(m_chartView, 1, 0);
    gridLayout->addWidget(m_selMapGroup, 0, 1);
    gridLayout->addWidget(m_mapGroup, 1, 1);

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addLayout(gridLayout);
    layout->addLayout(hLayout);
    setLayout(layout);

    setWindowTitle("Radiotherapy Simulator");
    setWindowIcon(QIcon("../Figures/logo.png"));

    QObject::connect(m_selMap, SIGNAL(currentIndexChanged(int)), this, SLOT(changeNumMap(int)));
    QObject::connect(m_selChart, SIGNAL(currentIndexChanged(int)), this, SLOT(changeChart(int)));
    QObject::connect(m_slider, SIGNAL(valueChanged(int)), m_timeS, SLOT(setValue(int)));
    QObject::connect(m_timeS, SIGNAL(valueChanged(int)), m_slider, SLOT(setValue(int)));
    QObject::connect(m_timeS, SIGNAL(valueChanged(int)), this, SLOT(changeIter(int)));
    QObject::connect(this, SIGNAL(updateSlider(int)), m_slider, SLOT(setValue(int)));
    QObject::connect(m_play, SIGNAL(clicked()), this, SLOT(play()));
    QObject::connect(m_sel, SIGNAL(clicked()), this, SLOT(sel()));
    QObject::connect(m_change, SIGNAL(clicked()), this, SLOT(change()));
    QObject::connect(m_newSim, SIGNAL(clicked()),this, SLOT(newSim()));
    QObject::connect(m_saveOutFiles, SIGNAL(clicked()), this, SLOT(saveOutFiles()));
    QObject::connect(m_saveChart, SIGNAL(clicked()), this, SLOT(saveChart()));
    QObject::connect(m_saveMap, SIGNAL(clicked()), this, SLOT(saveMap()));
    QObject::connect(m_saveAllMaps, SIGNAL(clicked()), this, SLOT(saveAllMaps()));

    drawMap(m_selMap->currentIndex(), m_slider->value());
    showMaximized();
}


void OutWindowOxy::change(){
    new InWindow;
    close();
}


void OutWindowOxy::changeChart(int numChart){
    switch(numChart){
    case 0:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_cHypDens->addSeries(m_sDash);
        m_sDash->attachAxis(m_xHypDens);
        m_sDash->attachAxis(m_yHypDens);
        m_chartView->setChart(m_cHypDens);
        break;
    }
    case 1:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_cPO2Stat->addSeries(m_sDash);
        m_sDash->attachAxis(m_xPO2Stat);
        m_sDash->attachAxis(m_yPO2Stat);
        m_chartView->setChart(m_cPO2Stat);
        m_chartView->chart()->legend()->markers(m_sDash).front()->setVisible(false);
        break;
    }

    case 2:{
        m_chartView->chart()->removeSeries(m_sDash);
        m_cVegfStat->addSeries(m_sDash);
        m_sDash->attachAxis(m_xVegfStat);
        m_sDash->attachAxis(m_yVegfStat);
        m_chartView->setChart(m_cVegfStat);
        m_chartView->chart()->legend()->markers(m_sDash).front()->setVisible(false);
        break;
    }
    }
    m_chartView->setRenderHint(QPainter::Antialiasing);
}


void OutWindowOxy::changeIter(int iter){
    drawMap(m_selMap->currentIndex(), iter);
    drawChartDashLine(iter);
}


void OutWindowOxy::changeNumMap(int numMap){
    drawMap(numMap, m_slider->value());
}


void OutWindowOxy::drawChartDashLine(int iter){
    QVector<QPointF> points;
    points.push_back(QPointF(iter, 0.0));
    points.push_back(QPointF(iter, 100.0));
    m_sDash->replace(points);
}


void OutWindowOxy::drawMap(int numMap, int mapIter){
    switch(numMap){
    case 0:{
        int iSclFac, jSclFac;
        for(int i(0); i < m_mapNrow; i++){
            iSclFac = i * m_mapSclFac;
            for(int j(0); j < m_mapNcol; j++){
                jSclFac = j * m_mapSclFac;
                for(int si(0); si < m_mapSclFac; si++){
                    for(int sj(0); sj < m_mapSclFac; sj++){
                        m_imPO2->setPixelColor(jSclFac + sj, iSclFac + si,
                                               QColor::fromHsv(240.0 * (m_maxpO2 - m_pO2[mapIter / m_simOxyTimeStep]
                                                               [i * m_mapNcol + j]) / m_maxpO2, 200, 255));
                    }
                }
            }
        }

        m_pixPO2->convertFromImage(*m_imPO2);
        m_map->setPixmap(m_pixPO2->scaledToWidth(1000));
        m_legendPO2->show();
        m_legendVegf->hide();
        break;
    }

    case 1:{
        int iSclFac, jSclFac;
        for(int i(0); i < m_mapNrow; i++){
            iSclFac = i * m_mapSclFac;
            for(int j(0); j < m_mapNcol; j++){
                jSclFac = j * m_mapSclFac;
                for(int si(0); si < m_mapSclFac; si++){
                    for(int sj(0); sj < m_mapSclFac; sj++){
                        m_imVegf->setPixelColor(jSclFac + sj, iSclFac + si,
                                                QColor::fromHsv(240.0 * (m_maxvegf - m_vegf[mapIter / m_simOxyTimeStep]
                                                                [i * m_mapNcol + j]) / m_maxvegf, 200, 255));
                    }
                }
            }
        }

        m_pixVegf->convertFromImage(*m_imVegf);
        m_map->setPixmap(m_pixVegf->scaledToWidth(1000));
        m_legendPO2->hide();
        m_legendVegf->show();
        break;
    }
    }

}


void OutWindowOxy::newSim(){
    new InWindow;
    close();
}


void OutWindowOxy::play(){
    m_stPlay = !m_stPlay;
    m_play->setIcon(QIcon("../Figures/pause.png"));
    for(int k(m_timeS->value()); k < m_timeS->maximum(); k++){
        if(!m_stPlay){
            m_play->setIcon(QIcon("../Figures/play.png"));
            return;
        }
        drawMap(m_selMap->currentIndex(), k);
        emit updateSlider(k);
        qApp->processEvents();
    }
    m_stPlay = !m_stPlay;
    m_play->setIcon(QIcon("../Figures/play.png"));
}


void OutWindowOxy::saveAllMaps(){
    QString fileName;
    QDir().mkdir("../Figures/pO2ms");
    const int K(m_pO2.size() * m_simOxyTimeStep);

    for(int k(0); k < K; k += m_simOxyTimeStep){
        drawMap(m_selMap->currentIndex(), k);
        emit updateSlider(k);
        fileName = QString::number(k);
        fileName = fileName.rightJustified(5, '0');
        fileName = "../Figures/pO2ms/" + fileName + ".png";
        m_map->grab().save(fileName);
        qApp->processEvents();
    }
}


void OutWindowOxy::saveChart(){
    QString fileName = QFileDialog::getSaveFileName(this, "Savem_chart", "../Figures",
                                                    "Images(*.png *.gif *.jpg *.jpeg)");
    if(!fileName.isEmpty()){
        m_chartView->chart()->series().back()->hide();
        m_chartView->grab().save(fileName);
        m_chartView->chart()->series().back()->show();
        m_chartView->chart()->legend()->markers().back()->setVisible(false);
    }
}


void OutWindowOxy::saveMap(){
    QString fileName = QFileDialog::getSaveFileName(this, "Save map", "../Figures",
                                                    "Images(*.png *.gif *.jpg *.jpeg)");
    if(!fileName.isEmpty()){
        m_map->grab().save(fileName);
    }
}


void OutWindowOxy::saveOutFiles(){
    QString dirName = QFileDialog::getSaveFileName(this, "Save ouptut files", "..");

    if(!dirName.isEmpty()){
        QDir newDir;
        newDir.mkdir(dirName);

        QDir outDir("../OutputFiles");
        QStringList filesList;
        filesList = outDir.entryList();

        for(int i(0); i < filesList.size(); i++){
            QFile::copy("../OutputFiles/" + filesList[i], dirName + "/" + filesList[i]);
        }
    }
}


void OutWindowOxy::sel(){
    new StartWindow;
    close();
}
