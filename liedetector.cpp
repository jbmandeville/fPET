#include <QDebug>
#include "liedetector.h"

lieDetector::lieDetector(const int numberSamples, const dVector BP0Values, const simEngine simulator, const PETRTM pet )
{
    _numberSamples  = numberSamples;
    _BP0Values      = BP0Values;
    _simulator      = simulator;
    _PETRTM         = pet;
}

void lieDetector::run()
{
    qDebug() << "lieDetector::run enter";
    double duration = _simulator.getDuration();
    double stepSize = _simulator.getStepSize();
    int lDownSample = _simulator.getDownSampling();
    int nTime = static_cast<int>(duration/stepSize) / lDownSample;
    dMatrix refRegion;      refRegion.resize(1);    refRegion[0].resize(nTime);
    dMatrix tissueVector;   tissueVector.resize(1); tissueVector[0].resize(nTime);

    _PETRTM.updateConditions();
    int nConditions = _PETRTM.getNumberConditions();
    if ( nConditions == 2 )
        _PETRTM.setCurrentCondition(1);

    int nSamplesTotal = _BP0Values.size() * _numberSamples;
    qDebug() << "lieDetector::run 3" << nSamplesTotal;
    int iDiv = qMax(1,nSamplesTotal/10);
    qDebug() << "lieDetector::run 3.5" << iDiv << _BP0Values.size() << _numberSamples;

    dMatrix errBPnd, errChall, tau2Ref;
    int nBPValues = _BP0Values.size();
    errBPnd.resize(nBPValues);   errChall.resize(nBPValues);  tau2Ref.resize(nBPValues);
    for (int jBP=0; jBP<nBPValues; jBP++)
    {
        errBPnd[jBP].resize(_numberSamples);    errChall[jBP].resize(_numberSamples);    tau2Ref[jBP].resize(_numberSamples);
        double BP0 = _BP0Values[jBP];
        _simulator.setBP0(BP0);
        for (int jSample=0; jSample<_numberSamples; jSample++)
        {
            int iSamplesTotal = jBP * _numberSamples + jSample;
            if ( iSamplesTotal%iDiv == 0 ) emit progressLieDetector(0);
            _simulator.run();  // run the simulations with randomized noise
            // Perform the analysis
            for (int jt=0; jt<nTime; jt++)
            {
                refRegion[0][jt]    = _simulator.getCrDown(jt);
                tissueVector[0][jt] = _simulator.getCtDown(jt);
            }
            _PETRTM.setReferenceRegion(refRegion);
            _PETRTM.setTissueVector(true,tissueVector);
            _PETRTM.prepare();
            _PETRTM.fitData(tissueVector);

            // update the BP error
            double guess = _PETRTM.getBP0InRun(0);
            errBPnd[jBP][jSample] = percentageError(guess,BP0);

            // update the challenge error
            double truth = _simulator.getChallengeMag();
            guess = getChallengeMagFromAnalysis();
            errChall[jBP][jSample] = guess - truth;

            // update the tau2Ref value
            if ( !_PETRTM.isRTM2() )
                guess = _PETRTM.getTau2RefInRun(0);
            else
                guess = 0.;
            tau2Ref[jBP][jSample] = guess;
        }
    }
    emit finishedLieDetector(errBPnd, errChall, tau2Ref);
    qDebug() << "lieDetector::run exit";
}

double lieDetector::getChallengeMagFromAnalysis()
{
    int nConditions = _PETRTM.getNumberConditions();
    if ( nConditions == 2 )
    {
        _PETRTM.evaluateCurrentCondition();
        double percentChange = _PETRTM.getBPndInCurrentCondition().x;
        percentChange *= 100./_PETRTM.getBP0InRun(0);
        return percentChange;
    }
    else
        return 0.;
}
