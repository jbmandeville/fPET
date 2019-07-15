#ifndef LIEDETECTOR_H
#define LIEDETECTOR_H

#include <QObject>
#include <QRunnable>
#include "simengine.h"
#include "petrtm.h"

class lieDetector: public QObject, public QRunnable
{
    Q_OBJECT   // this can create a problem of "missing vtable"; fix by "run qmake" from Qt Creator build menu
public:
    lieDetector(const int numberSamples, const dVector BP0Values, const simEngine simulator, const PETRTM pet );
    void run(); // if public, it can be run directly in the same thread to test thread safety
private:
    dVector _BP0Values;
    int _numberSamples;
    simEngine _simulator;
    PETRTM _PETRTM;

    inline double percentageError(double guess, double truth) {return 100.*(guess/truth-1.);}
    double getChallengeMagFromAnalysis();

signals:
    void progressLieDetector(int iProgress);
    void finishedLieDetector(dMatrix errBPnd, dMatrix errChall, dMatrix tau2Ref);
};

#endif // LIEDETECTOR_H
