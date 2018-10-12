#ifndef DEF_SIMTHREAD
#define DEF_SIMTHREAD

#include <QThread>

#include "../constOxyTissue.hpp"
#include "../coupler.hpp"
#include "../oxyTissue.hpp"
#include "../rootSimulator.hpp"
#include "../tissue.hpp"
#include "../treatment.hpp"

class SimThread : public QThread{
    Q_OBJECT
public:
    SimThread(QObject *parent);
    virtual void run();

signals:
    void progress(int value);
    void progressMax(int valueMax);
    void resultReady(int simType);
};

#endif
