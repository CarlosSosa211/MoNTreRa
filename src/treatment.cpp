/**
 * @file treatment.cpp
 * @brief
 * @author Carlos Sosa Marrero
 * @date 05.19.17
 */

#include <iostream>

#include "treatment.hpp"

using namespace std;

/*------------------------------------------------------------------------------
 * Constructor of the class Treatment.
 *
 * Inputs:
 *  - fraction: dose per fraction (Gy),
 *  - totalDose: total dose (Gy),
 *  - interval: interval between two sessions (h),
 *  - schedule: irradiation schedule (0, Mon-Fri; 1, Mon-Sund).
------------------------------------------------------------------------------*/

Treatment::Treatment(const double fraction, const double totalDose,
                     const double interval, const int schedule){
    int accSession(0), i(0), numSession;

    m_fraction  = fraction;
    m_interval  = interval;
    m_totalDose = totalDose;

    numSession = (int)(m_totalDose / m_fraction);

    switch(schedule){
    case 0:
        m_scheDescrip = "Mon-Fri";
        while(accSession <= numSession){
            if((i + 1) % 7 == 6 || (i + 1) % 7 == 0){
                m_schedule.push_back(false);
            }
            else{
                m_schedule.push_back(true);
                accSession++;
            }
            i++;
        }
        break;

    case 1:
        m_scheDescrip = "Mon-Sun";
        while(accSession <= numSession){
            m_schedule.push_back(true);
            accSession++;
        }
        break;

    default:
        m_scheDescrip = "Mon-Fri";
        while(accSession <= numSession){
            if((i + 1) % 7 == 6 || (i + 1) % 7 == 0){
                m_schedule.push_back(false);
            }
            else{
                m_schedule.push_back(true);
                accSession++;
            }
            i++;
        }
        break;
    }
}


/*------------------------------------------------------------------------------
 * Destructor of the class Treatment.
------------------------------------------------------------------------------*/

Treatment::~Treatment(){
}


/*------------------------------------------------------------------------------
 * This function gets the duration of the treatment.
------------------------------------------------------------------------------*/

double Treatment::getDuration() const{
    return m_schedule.size() * m_interval;
}


/*------------------------------------------------------------------------------
 * This function gets the fraction.
------------------------------------------------------------------------------*/

double Treatment::getFraction() const{
    return m_fraction;
}


/*------------------------------------------------------------------------------
 * This function gets the interval between two sessions.
------------------------------------------------------------------------------*/

double Treatment::getInterval() const{
    return m_interval;
}


/*------------------------------------------------------------------------------
 * This function gets the schedule.
------------------------------------------------------------------------------*/

vector<bool> Treatment::getSchedule() const{
    return m_schedule;
}

/*------------------------------------------------------------------------------
 * This function gets the schedule description.
------------------------------------------------------------------------------*/

string Treatment::getScheDescrip() const{
    return m_scheDescrip;
}


/*------------------------------------------------------------------------------
 * This function gets the total dose.
------------------------------------------------------------------------------*/

double Treatment::getTotalDose() const{
    return m_totalDose;
}


/*------------------------------------------------------------------------------
 * Overload of the << operator.
------------------------------------------------------------------------------*/

ostream &operator<<(ostream &stream, Treatment *treatment){
    if(treatment){
        stream << "Total dose = " << treatment->getTotalDose() << " Gy" << endl;
        stream << "Fraction = " << treatment->getFraction() << " Gy" << endl;
        stream << "Interval = " << treatment->getInterval() << " h" << endl;
        stream << "Schedule = " << treatment->getScheDescrip() <<  endl;
        stream << "---------------------------------------------" << endl;
    }

    else{
        stream << "No treatment" <<  endl;
        stream << "---------------------------------------------" << endl;
    }
    return stream;
}
