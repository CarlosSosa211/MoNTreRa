HEADERS += \
    ../absOxyCell.hpp \
    ../cell.hpp \
    ../constOxyCell.hpp \
    ../constOxyTissue.hpp \
    ../coupler.hpp \
    inWindow.hpp \
    ../model.hpp \
    outWindow.hpp \
    outWindow3D.hpp \
    outWindowOxy.hpp \
    ../oxyCell.hpp \
    ../oxyTissue.hpp \
    ../rootSimulator.hpp \
    ../simpcell.hpp \
    ../simulator.hpp \
    startWindow.hpp \
    ../tissue.hpp \
    ../treatment.hpp \
    simThread.hpp \
    inAnWindow.hpp \
    style.hpp \
    switch.hpp \
    outAnWindow.hpp


QT += \
    charts \
    #datavisualization \
    widgets


SOURCES += \
    ../absOxyCell.cpp \
    ../cell.cpp \
    ../constOxyCell.cpp \
    ../constOxyTissue.cpp \
    ../coupler.cpp \
    inWindow.cpp \
    main.cpp \
    ../model.cpp \
    outWindow.cpp \
    outWindow3D.cpp \
    outWindowOxy.cpp \
    ../oxyCell.cpp \
    ../oxyTissue.cpp \
    ../rootSimulator.cpp \
    ../simulator.cpp \
    startWindow.cpp \
    ../tissue.cpp \
    ../treatment.cpp \
    simThread.cpp \
    inAnWindow.cpp \
    switch.cpp \
    outAnWindow.cpp

DISTFILES += \
    ../InputFiles/in.dat \
    ../InputFiles/inAng.dat \
    ../InputFiles/inConstOxy.dat \
    ../InputFiles/inOxy.dat \
    ../InputFiles/inTG.dat \
    ../InputFiles/inSlowSim.dat \
    ../InputFiles/inTestOxy.dat
