#ifndef LIEDETECTOR_H
#define LIEDETECTOR_H

#include <QObject>
#include <QRunnable>
#include "simengine.h"
#include "petrtm.h"

class lieDetectorBPnd: public QObject, public QRunnable
{
    Q_OBJECT   // this can create a problem of "missing vtable"; fix by "run qmake" from Qt Creator build menu
public:
    lieDetectorBPnd(const int numberSamples, const dVector BP0Values, const simEngine simulator, const PETRTM pet );
    void run(); // if public, it can be run directly in the same thread to test thread safety
private:
    dVector _BP0Values;
    int _numberSamples;
    simEngine _simulator;
    bool _fitk4;
    PETRTM _PETRTM;

    inline double percentageError(double guess, double truth) {return 100.*(guess/truth-1.);}
    double getChallengeMagFromAnalysis();

signals:
    void progresslieDetector(int iProgress);
    void finishedlieDetector(dMatrix errBPnd, dMatrix errChall, dMatrix tau2Ref, dMatrix errTau4, dMatrix errDV, double sigma2);
};

class lieDetectorTau4: public QObject, public QRunnable
{
    Q_OBJECT   // this can create a problem of "missing vtable"; fix by "run qmake" from Qt Creator build menu
public:
    lieDetectorTau4(const int numberSamples, const dVector tau4Values, const simEngine simulator, const PETRTM pet );
    void run(); // if public, it can be run directly in the same thread to test thread safety
private:
    dVector _tau4Values;
    int _numberSamples;
    simEngine _simulator;
    bool _fitk4;
    PETRTM _PETRTM;

    inline double percentageError(double guess, double truth) {return 100.*(guess/truth-1.);}
    double getChallengeMagFromAnalysis();

signals:
    void progresslieDetector(int iProgress);
    void finishedlieDetector(dMatrix errBPnd, dMatrix errChall, dMatrix AIC, double sigma2);
};

#endif // LIEDETECTOR_H
