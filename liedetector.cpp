#include <QDebug>
#include "liedetector.h"

lieDetectorBPnd::lieDetectorBPnd(const int numberSamples, const dVector BP0Values, const simEngine simulator, const PETRTM pet )
{
    _numberSamples  = numberSamples;
    _BP0Values      = BP0Values;
    _simulator      = simulator;
    _PETRTM         = pet;
    _fitk4          = _PETRTM.getFitk4() && _PETRTM.isForwardModel();
}

void lieDetectorBPnd::run()
{
    FUNC_ENTER;
    int nTime = _simulator.getNumberBins();
    dMatrix tissueVector;   tissueVector.resize(1); tissueVector[0].resize(nTime);

    _PETRTM.updateConditions();
    int nConditions = _PETRTM.getNumberConditions();
    if ( nConditions == 2 )
        _PETRTM.setCurrentCondition(1);

    int nBPValues = _BP0Values.size();
    int nSamplesTotal = nBPValues * _numberSamples;
    int iDiv = qMax(1,nSamplesTotal/10);

    double tau4Guess = _PETRTM._simulator[0].getTau4Nominal();
    if ( tau4Guess == 0. ) tau4Guess = 10.;

    dMatrix errBPnd, errChall, tau2Ref, errTau4;
    errBPnd.resize(nBPValues);   errChall.resize(nBPValues);  tau2Ref.resize(nBPValues);  errTau4.resize(nBPValues);
    double sigma2Sum = 0.;
    for (int jBP=0; jBP<nBPValues; jBP++)
    {
        errBPnd[jBP].resize(_numberSamples); errChall[jBP].resize(_numberSamples); tau2Ref[jBP].resize(_numberSamples); errTau4[jBP].resize(_numberSamples);
        double BP0 = _BP0Values[jBP];
        _simulator.setBP0(BP0);
        for (int jSample=0; jSample<_numberSamples; jSample++)
        {
            int iSamplesTotal = jBP * _numberSamples + jSample;
            if ( iSamplesTotal%iDiv == 0 ) emit progresslieDetector(0);
            _simulator.run();  // run the simulations with randomized noise
            // Perform the analysis
            for (int jt=0; jt<nTime; jt++)
                tissueVector[0][jt] = _simulator.getCtCoarse(jt);
            _PETRTM.setTissueRegion(tissueVector);
            _PETRTM.prepare();
            _PETRTM.fitData(tissueVector);

            if ( _PETRTM.isForwardModel() )
                sigma2Sum += _PETRTM.getSimulationSigma2(0);
            else
                sigma2Sum += _PETRTM.getSigma2();

            // update the BP error
            double guess = _PETRTM.getBP0InRun(0).x;
            errBPnd[jBP][jSample] = percentageError(guess,BP0);

            // update the challenge error
            double truth = _simulator.getChallengeMag();
            guess = getChallengeMagFromAnalysis();
            errChall[jBP][jSample] = guess - truth;

            // update the k4 error
            truth = _simulator.getTau4();
            guess = _PETRTM.getTau4InRun(0);
            errTau4[jBP][jSample] = percentageError(guess,truth);

            // update the tau2Ref value
            if ( !_PETRTM.isRTM2() )
                guess = _PETRTM.getTau2RefInRun(0);
            else
                guess = 0.;
            tau2Ref[jBP][jSample] = guess;
        }
    }
    sigma2Sum /= static_cast<double>(nBPValues*_numberSamples);
    emit finishedlieDetector(errBPnd, errChall, tau2Ref, errTau4, sigma2Sum);
    FUNC_EXIT;
}

double lieDetectorBPnd::getChallengeMagFromAnalysis()
{
    int nConditions = _PETRTM.getNumberConditions();
    if ( nConditions == 2 )
    {
        _PETRTM.evaluateCurrentCondition();
        double percentChange = _PETRTM.getBPndInCurrentCondition();
        percentChange *= 100./_PETRTM.getBP0InRun(0).x;
        return percentChange;
    }
    else
        return 0.;
}

lieDetectorTau4::lieDetectorTau4(const int numberSamples, const dVector tau4Values, const simEngine simulator, const PETRTM pet )
{
    _numberSamples  = numberSamples;
    _tau4Values      = tau4Values;
    _simulator      = simulator;
    _PETRTM         = pet;
    _fitk4          = _PETRTM.getFitk4() && _PETRTM.isForwardModel();
}

void lieDetectorTau4::run()
{
    FUNC_ENTER;
    int nTime = _simulator.getNumberBins();
    dMatrix tissueVector;   tissueVector.resize(1); tissueVector[0].resize(nTime);

    _PETRTM.updateConditions();
    int nConditions = _PETRTM.getNumberConditions();
    if ( nConditions == 2 )
        _PETRTM.setCurrentCondition(1);

    int nTau4Values = _tau4Values.size();
    int nSamplesTotal = nTau4Values * _numberSamples;
    int iDiv = qMax(1,nSamplesTotal/10);

    double tau4Guess = _PETRTM._simulator[0].getTau4Nominal();
    if ( tau4Guess == 0. ) tau4Guess = 10.;

    dMatrix errBPnd, errChall, AIC;
    errBPnd.resize(nTau4Values);   errChall.resize(nTau4Values);  AIC.resize(nTau4Values);
    double sigma2Sum = 0.;
    for (int jBP=0; jBP<nTau4Values; jBP++)
    {
        errBPnd[jBP].resize(_numberSamples); errChall[jBP].resize(_numberSamples);  AIC[jBP].resize(_numberSamples);
        double tau4 = _tau4Values[jBP];
        _PETRTM.setTau4(0,tau4);
        for (int jSample=0; jSample<_numberSamples; jSample++)
        {
            int iSamplesTotal = jBP * _numberSamples + jSample;
            if ( iSamplesTotal%iDiv == 0 ) emit progresslieDetector(0);
            _simulator.run();  // run the simulations with randomized noise
            // Perform the analysis
            for (int jt=0; jt<nTime; jt++)
                tissueVector[0][jt] = _simulator.getCtCoarse(jt);
            _PETRTM.setTissueRegion(tissueVector);
            _PETRTM.prepare();
            _PETRTM.fitData(tissueVector);

            sigma2Sum += _PETRTM.getSimulationSigma2(0);

            // update the BP error
            double truth = _simulator.getBP0();
            double guess = _PETRTM.getBP0InRun(0).x;
            errBPnd[jBP][jSample] = percentageError(guess,truth);

            // update the challenge error
            truth = _simulator.getChallengeMag();
            guess = getChallengeMagFromAnalysis();
            errChall[jBP][jSample] = guess - truth;

            // update the k4 error
            truth = _simulator.getTau4();
            guess = _PETRTM.getTau4InRun(0);
            AIC[jBP][jSample] = percentageError(guess,truth);
        }
    }
    sigma2Sum /= static_cast<double>(nTau4Values*_numberSamples);
    emit finishedlieDetector(errBPnd, errChall, AIC, sigma2Sum);
    FUNC_EXIT;
}

double lieDetectorTau4::getChallengeMagFromAnalysis()
{
    int nConditions = _PETRTM.getNumberConditions();
    if ( nConditions == 2 )
    {
        _PETRTM.evaluateCurrentCondition();
        double percentChange = _PETRTM.getBPndInCurrentCondition();
        percentChange *= 100./_PETRTM.getBP0InRun(0).x;
        return percentChange;
    }
    else
        return 0.;
}
