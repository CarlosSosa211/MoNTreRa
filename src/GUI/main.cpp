/**
 * @file main.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @author Nicolas Ciferri
 * @author Alfredo Hernandez
 * @date 05.19.17
 */

#include <QApplication>

#include "inWindow.hpp"
#include "outWindow.hpp"
#include "outWindow3D.hpp"
#include "startWindow.hpp"

using namespace std;

int main(int argc, char *argv[]){
    QApplication app(argc, argv);
    StartWindow startWindow;
    return app.exec();
}
