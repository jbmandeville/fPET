#include "simanalyzer.h"
#include <QtMath>
#include <QDebug>
//#include <QRandomGenerator> // include only in Qt5.10 and higher

#define RAND_0_1 (rand()/2147483647.)

simAnalyzer::simAnalyzer()
{
    // initialize
    initializeBins();
    // run
    _parValues = {_BP0,_R1,_k2Ref,_tau4,_dBPndChallenge, _plasmaMag, _dBPndUptake};
    if ( ! _RTM3 ) _parValues[2] = _k2RefFixed;
    runAndCalculateCost();
}
void simAnalyzer::initializeBins()
{
    static int NSAMPLESPERBIN = 10;
    if ( _nBins != _durationBinSec.size() )
    {
        double defaultBinSize=60;
        if ( _durationBinSec.size() > 0 ) defaultBinSize = _durationBinSec[0];
        _durationBinSec.fill(60,_nBins);  // duration in min
    }
    if ( _nBins != _numberSamplesPerBin.size() )
        _numberSamplesPerBin.fill(NSAMPLESPERBIN,_nBins);
    updateFineSamples();
}
void simAnalyzer::updateFineSamples()
{
    _dtFine.clear();  _timeFine.clear(); _timeCoarse.clear();  _dtCoarse.clear();
    double time=0.;
    for (int jBin=0; jBin<_nBins; jBin++)
    {
        double dt = static_cast<double>(_durationBinSec[jBin])/60.;
        double dtFine = dt / _numberSamplesPerBin[jBin];
        double avTime=0.;
        for (int jSample=0; jSample<_numberSamplesPerBin[jBin]; jSample++)
        {
            _dtFine.append(dtFine);
            _timeFine.append(time + dtFine/2.);
            avTime += _timeFine.last();
            time += dtFine;  // start of bin
        }
        avTime /= static_cast<double>(_numberSamplesPerBin[jBin]);
        _timeCoarse.append(avTime);
        _dtCoarse.append(dt);
    }
    if ( _weights.size() != _dtCoarse.size() ) _weights.fill(1.,_dtCoarse.size());
    updateChallengeShape();
}

void simAnalyzer::setCrFit(dVector CrFitCoarse)
{
    setNumberBins(CrFitCoarse.size());
    _CrBinned = CrFitCoarse;

    tk::spline spline;
    spline.set_points(_timeCoarse,_CrBinned);    // currently it is required that X is already sorted

    // fine scale: interpolate LOESS using spline
    int nTimeFine = getNumberTimeBinsFine();
    _Cr.clear();
    for (int jTime=0; jTime<nTimeFine; jTime++)
        _Cr.append(qMax(0.,spline(_timeFine[jTime])));

    _CrDot.fill(0.,_Cr.size());
    // Calculate the derivative
    for (int jBin=0; jBin<_Cr.size(); jBin++)
    {
        if ( jBin != 0 )
            _CrDot[jBin] = (_Cr[jBin] - _Cr[jBin-1]) / _dtFine[jBin];
        else
            _CrDot[jBin] = 0.;
    }
}

double simAnalyzer::runAndCalculateCost()
{
//    FUNC_ENTER << _parValues;
    if ( _CtData.size() == 0. ) return 0.;

    _BP0   = _parValues.at(0);
    _R1    = _parValues.at(1);
    _k2Ref = _parValues.at(2);
    _tau4  = _parValues.at(3);
    _dBPndChallenge = _parValues.at(4);
    _plasmaMag      = _parValues.at(5);
    _dBPndUptake    = _parValues.at(6);
    if ( _tau4 < 1. )
        generateTargetTACSRTM();
    else
        generateTargetTACFRTM();
    _sigma2 = 0.;
    int nTime = _CtData.size();
    for (int jt=0; jt<nTime; jt++)
    {
//        if ( _weight_t[jt] != 0. )
//        {
//        FUNC_INFO << "weight[" << jt << "] =" << _weights[jt];
        _sigma2 += _weights[jt] * SQR(_CtData[jt]-_CtBinned[jt]);
//        }
    }
    double dof = static_cast<double>(nTime);
    _sigma2 /= dof;
//    FUNC_EXIT << "dof" << dof << "sigma2" << _sigma2;
    if ( isnan(_sigma2) || isinf(_sigma2) )
        return 1.e300;
    else
        return _sigma2;
}

void simAnalyzer::generateTargetTACFRTM()
{
    // dCr_dt = K1Ref * Cp - k2Ref * Cr
    // dCb_dt = k3 * Cf - koff * Cb = kon * Bavail * Cf - koff * Cb
    // dCf_dt = K1 * Cp + koff * Cb - (k2+k3)*Cf
    // R1 = K1/K1Ref = k2/k2Ref

    if ( _dtFine.size() == 0 ) return;
    if ( _Cr.size()     == 0 ) return;
    if ( _CrDot.size()  == 0 ) return;

    // Derived quantities: update tissue properties
    _k2  = _R1 * _k2Ref;
    if ( _tau4 == 0. ) qFatal("Error: tau4 == 0 in generateTargetTACFRTM(); use SRTM.");
    double k4 = 1./_tau4;
    _k3  = _BP0  * k4;

    double k3=_k3;
    double Cf=0.;  double Cb=0.;  double lastCr=0.;
    int nTime = _dtFine.size();

    _Ct.clear();
    //        double Cp = 0;        // plasma contribution
    for ( int jt=0; jt<nTime; jt++)
    {
        int iUptakeTime       = static_cast<int>(_uptakeTime    / _dtFine[jt]);
        bool uptake = (jt < iUptakeTime);
        if ( uptake && _dBPndUptake != 0. )
        {
            double BP1 = _BP0 - _dBPndUptake;
            k3 = BP1 * k4;
        }
        else
        {
            double BP1 = _BP0 - _dBPndChallenge * _challengeCurve[jt];
            k3 = BP1 * k4;
        }
        double Cr = _Cr[jt];
        // To add plasma contamination, one needs K1Rref (not just R1)
        double Cp = 0.;
        if ( jt > 0 )
        {
            Cp = _plasmaMag * ( _k2Ref * _Cr[jt] + _CrDot[jt]);        // plasma
            double dCfdt = _R1 * ( _k2Ref * _Cr[jt] + _CrDot[jt]) + k4 * Cb - (_k2 + k3) * Cf;
            double dCbdt = k3 * Cf - k4 * Cb;
            Cf += dCfdt * _dtFine[jt];
            Cb += dCbdt * _dtFine[jt];
            lastCr = Cr;
        }
        _Ct.append(Cf + Cb + Cp);
        lastCr = Cr;
    }
    // Downsample the TAC
    _CtBinned = downSample(_Ct);
}

void simAnalyzer::generateTargetTACSRTM()
{
    // dCr_dt = K1Ref * Cp - k2Ref * Cr  => Cp = (dCr_dt + k2Ref * Cr)/K1Ref
    // dCt_dt = K1 * Cp - k2a * Ct = R1 * (dCr_dt + k2Ref * Cr) - k2a * Cr
    // R1 = K1/K1Ref = k2/k2Ref

    // Derived quantities: update tissue properties
    _k2  = _R1 * _k2Ref;
    double k2a = _k2 / (1. + _BP0);

    //        double Cp = 0;        // plasma contribution
    double lastCr=0.;
    int nTime = _dtFine.size();
    _Ct.clear();  double Ct = 0.;
    for ( int jt=0; jt<nTime; jt++)
    {
        int iUptakeTime = static_cast<int>(_uptakeTime / _dtFine[jt]);
        bool uptake = (jt < iUptakeTime);
        if ( uptake && _dBPndUptake != 0. )
        {
            double BP1 = _BP0 - _dBPndUptake;
            k2a = _k2 / (1. + BP1);
        }
        else
        {
            double BP1 = _BP0 - _dBPndChallenge * _challengeCurve[jt];
            k2a = _k2 / (1. + BP1);
        }
        double Cr = _Cr[jt];
        // To add plasma contamination, one needs K1Rref (not just R1)
        double Cp = 0.;
        if ( jt > 0 )
        {
            Cp = _plasmaMag * ( _k2Ref * _Cr[jt] + _CrDot[jt]);
            double dCtdt = _R1 * ( _k2Ref * _Cr[jt] + _CrDot[jt]) - k2a * Ct;
            Ct += dCtdt * _dtFine[jt];
            lastCr = Cr;
        }
        _Ct.append(Ct + Cp);
        lastCr = Cr;
    }
    // Downsample the TAC
    _CtBinned = downSample(_Ct);
}

dVector simAnalyzer::downSample(dVector original)
{
    if ( original.size() != _dtFine.size() )
        qFatal("Error: the fine-scale TAC does not match the fine-scale bin length in downSample().");
    dVector binned;
    int iTFine = 0;
    for (int jBin=0; jBin<_nBins; jBin++)
    {
        double binValue = 0.;
        for (int jt=0; jt<_numberSamplesPerBin[jBin]; jt++, iTFine++)
            binValue += original[iTFine];
        binned.append(binValue/ static_cast<double>(_numberSamplesPerBin[jBin]));
    }
    return binned;
}

double simAnalyzer::getDurationScan()
{
    double duration=0.;
    for (int jt=0; jt<_dtFine.size(); jt++)
        duration += _dtFine[jt];
    return duration;
}
double simAnalyzer::getdk2a()
{
    // k2a = k2 / (1+_BP0)
    // BP1 = BP0
    double k2a_0 = _k2/(1.+_BP0);
    double BP1 = _BP0 - _dBPndChallenge;
    double k2a_1 = _k2/(1.+BP1);
    double dk2a = k2a_1 - k2a_0;
    return dk2a;
}
double simAnalyzer::getdk2k3()
{
    double BP1 = _BP0 - _dBPndChallenge;
    double DBP = _BP0 - BP1;
    return - _k2 * DBP;  // don't use PET convention
}

void simAnalyzer::fit(double BP0, double R1, double k2Ref, double dBPnd, double tau4,
                      bool fitk4, bool RTM3, bool challenge, int nBootStrap)
{
    FUNC_ENTER << tau4 << fitk4;
//    FUNC_ENTER << "BP0" << _BP0 << "R1" << _R1 << "k2Ref" << k2Ref << "dBPnd" << dBPnd << "bool RTM3" << RTM3 << "bool challenge" << challenge;
//    FUNC_ENTER << "BP5, R1, k2Ref = " << BP0 << R1 << k2Ref << dBPnd << RTM3 << challenge;
    _RTM3 = RTM3;
    _BP0   = BP0;
    _R1    = R1;
    _tau4  = tau4;
    if ( RTM3 )
        _k2Ref = k2Ref;
    else
        _k2Ref = _k2RefFixed;
    _fitk4 = fitk4;
    if ( fitk4 )
        _tau4 = tau4;
    _dBPndChallenge = dBPnd;
    _dBPndUptake = _plasmaMag = 0.;
    bool uptake = false;
    _parAdjust = {true, true, RTM3, _fitk4, challenge, _includePlasma, uptake};
    fitTACByGridSearch(6);  // 6 -> 1
    fitTACByLineScan(5);
    fitTACByLineScan(4);
    fitTACByLineScan(3);
    fitTACByLineScan(2);
    fitTACByLineScan(1);
    if ( nBootStrap > 0 )
        calculateErrorsByBootstrapping(nBootStrap);
    else
        _BP0Err = _R1Err = _k2RefErr = _dBPndErr = 0.;
//    FUNC_EXIT << "BP0, R1, k2Ref, dBPnd = " << _BP0 << _R1 << _k2Ref << _dBPndChallenge << _dBPndUptake;
}

void simAnalyzer::calculateErrorsByBootstrapping(int nBootStrap)
{
    // save original values
    double BP0Save   = getBP0();
    double R1Save    = getR1();
    double k2RefSave = getk2Ref();
    double dBPndSave = getChallengeMag();
    dVector fitSave  = _CtBinned;
    //
    int nTime = _CtData.size();
    dVector residue;  residue.fill(0.,nTime);
    double tau = 20.2334 * 1.442695;    // convert from half life to exponential time constant
    for (int jt=0; jt<nTime; jt++)
    {
        double dt = static_cast<double>(_dtCoarse[jt])/60.;
        double time = _timeCoarse[jt];
        double correction = qExp(time/tau) / dt;
        residue[jt] = (_CtData[jt]-fitSave[jt]) * correction;
    }

    // Perform a large number of trials using random replacement of residual errors
    dVector BP0Trial, R1Trial, k2RefTrial, dBPTrial;
    dVector CtDataSynthetic;  CtDataSynthetic.fill(0.,_CtData.size());
    for (int jTrial=0; jTrial<nBootStrap; jTrial++)
    {
        for (int jt=0; jt<nTime; jt++)
        {
            double dt = static_cast<double>(_dtCoarse[jt])/60.;
            double time = _timeCoarse[jt];
            double correction = qExp(time/tau) / dt;
            int iRand = nTime * RAND_0_1;
            CtDataSynthetic[jt] = fitSave[jt] + residue[iRand] / correction;
        }
        setCtData(CtDataSynthetic);
//        fit(BP0Save, R1Save, k2RefSave, dBPndSave, _RTM3, true);
        _BP0   = BP0Save;
        _R1    = R1Save;
        _k2Ref = k2RefSave;
        _dBPndChallenge = dBPndSave;
        fitTACByLineScan(3);
        fitTACByLineScan(2);
        fitTACByLineScan(1);
        BP0Trial.append(getBP0());
        R1Trial.append(getR1());
        k2RefTrial.append(getk2Ref());
        dBPTrial.append(getChallengeMag());
    }

    // compute averages
    double BP0Av = 0.;  double R1Av = 0.;  double k2RefAv = 0.;  double dBPAv = 0.;
    for (int jTrial=0; jTrial<nBootStrap; jTrial++)
    {
        BP0Av   += BP0Trial[jTrial];
        R1Av    += R1Trial[jTrial];
        k2RefAv += k2RefTrial[jTrial];
        dBPAv   += dBPTrial[jTrial];
    }
    double trials = nBootStrap;
    BP0Av /= trials;  R1Av /= trials;  k2RefAv /= trials;   dBPAv /= trials;

    // compute standard deviations
    double BP0SEM = 0.;  double R1SEM = 0.;  double k2RefSEM = 0.;  double dBPSEM = 0.;
    for (int jTrial=0; jTrial<nBootStrap; jTrial++)
    {
        BP0SEM   += SQR(BP0Av - BP0Trial[jTrial]);
        R1SEM    += SQR(R1Av - R1Trial[jTrial]);
        k2RefSEM += SQR(k2RefAv - k2RefTrial[jTrial]);
        dBPSEM   += SQR(dBPAv - dBPTrial[jTrial]);
    }
    BP0SEM /= trials-1.;      R1SEM /= trials-1.;    k2RefSEM /= trials-1.;       dBPSEM /= trials-1.;
    BP0SEM = qSqrt(BP0SEM);   R1SEM = qSqrt(R1SEM);  k2RefSEM = qSqrt(k2RefSEM);  dBPSEM = qSqrt(dBPSEM);
    _BP0Err = BP0SEM;
    _R1Err  = R1SEM;
    _k2RefErr = k2RefSEM;
    _dBPndErr = dBPSEM;

    FUNC_INFO << "BP0" << BP0Save << BP0Av << "±" << _BP0Err << "dBP" << dBPndSave <<  dBPAv << "±" << _dBPndErr;

    // restore original values
    _BP0   = BP0Save;
    _R1    = R1Save;
    _k2Ref = k2RefSave;
    _dBPndChallenge = dBPndSave;
    _CtBinned = fitSave;
}

void simAnalyzer::fitTACByGridSearch(int level)
{
//    FUNC_ENTER << "BP0, R1, k2Ref, dBPnd = " << _BP0 << _R1 << _k2Ref << _dBPndChallenge << _dBPndUptake << "level" << level;
    // Fit a gamma function plus a polynomial by repeated GLM using a grid search for the best gamma parameters (3)
    double BP0    = _BP0;
    double R1     = _R1;
    double k2Ref  = _k2Ref;
    double tau4   = _tau4;
    double dBPndC = _dBPndChallenge;
    double dBPndU = _dBPndUptake;
    double plasma = _plasmaMag;

//    double BP0Width = qMax(_BP0/5,1.);
//    double R1Width = qMax(_R1/2.,0.5);
    double BP0Width = qMax(_BP0,1.);
    double R1Width = qMax(_R1,0.5);
    double k2RefWidth = qMax(_k2Ref,.2);
    double tau4Width  = qMax(_tau4,10.);
    double dBPndCWidth = qMax(_BP0/2.,0.1);
    double dBPndUWidth = qMax(_BP0/2.,0.1);
    double plasmaWidth = qMax(_plasmaMag/2.,0.001);

    if ( level == 5 )
    {
        BP0Width /= 2.;
        k2RefWidth /= 2.;
        tau4Width /= 2.;
        R1Width /= 2.;
        dBPndCWidth /= 2.;
        dBPndUWidth /= 2.;
        plasmaWidth /= 2.;
    }
    else if ( level == 4 )
    {
        BP0Width /= 4.;
        k2RefWidth /= 4.;
        tau4Width /= 4.;
        R1Width /= 4.;
        dBPndCWidth /= 4.;
        dBPndUWidth /= 4.;
        plasmaWidth /= 4.;
    }
    else if ( level == 3 )
    {
        BP0Width/=8;
        k2RefWidth /= 8.;
        tau4Width /= 8.;
        R1Width /= 8.;
        dBPndCWidth /= 8.;
        dBPndUWidth /= 8.;
        plasmaWidth /= 8.;
    }
    else if ( level == 2 )
    {
        BP0Width/=16.;
        k2RefWidth /= 16.;
        tau4Width /= 16.;
        R1Width /= 16.;
        dBPndCWidth /= 16.;
        dBPndUWidth /= 16.;
        plasmaWidth /= 16.;
    }
    else // if ( level == 1 )
    {
        BP0Width/= 32.;
        k2RefWidth /= 32.;
        tau4Width /= 32.;
        R1Width /= 32.;
        dBPndCWidth /= 32.;  // 3%
        dBPndUWidth /= 32.;  // 3%
        plasmaWidth /= 32.;
    }

//    FUNC_INFO << "widths for BP0Width, R1Width, k2RefWidth, dBPndCWidth" << BP0Width << R1Width << k2RefWidth << dBPndCWidth;

    double stepDivider = 5.;
    double BP0Step   = BP0Width  /stepDivider;
    double R1Step    = R1Width   /stepDivider;
    double k2RefStep = k2RefWidth/stepDivider;
    double tau4Step  = tau4Width/stepDivider;
    double dBPndCStep = dBPndCWidth/stepDivider;
    double dBPndUStep = dBPndUWidth/stepDivider;
    double plasmaStep = plasmaWidth/stepDivider;

    double BP0Start   = qMax(0.,BP0-BP0Width);      double BP0Stop = BP0 + BP0Width;
    double R1Start    = qMax(0.,R1-R1Width);        double R1Stop   = R1 + R1Width;
    double k2RefStart = qMax(0.,k2Ref-k2RefWidth);  double k2RefStop = k2Ref + k2RefWidth;
    double tau4Start  = qMax(0.,tau4-tau4Width);    double tau4Stop = tau4 + tau4Width;
    double dBPndCStart = dBPndC-dBPndCWidth;        double dBPndCStop = dBPndC + dBPndCWidth;
    double dBPndUStart = dBPndU-dBPndUWidth;        double dBPndUStop = dBPndU + dBPndUWidth;
    double plasmaStart = plasma-plasmaWidth;        double plasmaStop = plasma + plasmaWidth;

    _parValues     = {_BP0, _R1, _k2Ref, _tau4, _dBPndChallenge, _plasmaMag, _dBPndUptake};
    _parIncrements = {BP0Step, R1Step, k2RefStep, tau4Step, dBPndCStep, plasmaStep, dBPndUStep};
    dVector bestSet = _parValues;
    double bestSigma2 = runAndCalculateCost();

    if ( !_parAdjust[0] )
        BP0Start = BP0Step = _BP0;
    if ( !_parAdjust[1] )
        R1Start = R1Stop = _R1;
    if ( !_parAdjust[2] )
        k2RefStart = k2RefStop = _k2Ref;
    if ( !_parAdjust[3] )
        tau4Start = tau4Stop = _tau4;
    if ( !_parAdjust[4] )
        dBPndCStart = dBPndCStop = _dBPndChallenge;
    if ( !_parAdjust[5] )
        plasmaStart = plasmaStop = _plasmaMag;
    if ( !_parAdjust[6] )
        dBPndUStart = dBPndUStop = _dBPndUptake;

    for ( BP0=BP0Start; BP0<BP0Stop; BP0+=BP0Step)
    {
        _BP0 = BP0;
        for ( R1=R1Start; R1<R1Stop; R1+=R1Step)
        {
            _R1 = R1;
            for ( k2Ref=k2RefStart; k2Ref<k2RefStop; k2Ref+=k2RefStep)
            {
                _k2Ref = k2Ref;
                for ( tau4=tau4Start; tau4<tau4Stop; tau4+=tau4Step)
                {
                    _tau4 = tau4;
                    for ( dBPndC=dBPndCStart; dBPndC<dBPndCStop; dBPndC+=dBPndCStep)
                    {
                        _dBPndChallenge = dBPndC;
                        for ( plasma=plasmaStart; plasma<plasmaStop; plasma+=plasmaStep)
                        {
                            _plasmaMag = plasma;
                            for ( dBPndU=dBPndUStart; dBPndU<dBPndUStop; dBPndU+=dBPndUStep)
                            {
                                _dBPndUptake = dBPndU;
                                _parValues = {_BP0,_R1,_k2Ref,_tau4,_dBPndChallenge, _plasmaMag, _dBPndUptake};
                                double sigma2 = runAndCalculateCost();
                                if ( sigma2 < bestSigma2 )
                                {
                                    bestSigma2 = sigma2;
                                    bestSet = _parValues;
                                }
                            }  // dBPndU
                        } // plasmaMag
                    }  // dBPndC
                } // tau4
            } // k2Ref
        } // R1
    } // BP0
    _parValues = bestSet;
    bestSigma2 = runAndCalculateCost();
//    FUNC_EXIT << "bestSet " << bestSet << "with sigma =" << qSqrt(bestSigma2) << "and tau4 = " << _tau4 << "\n";
}

void simAnalyzer::fitTACByLineScan(int level) // toleranceCost is a fraction (e.g., 0.01 is 1%)
{
//    FUNC_ENTER << "BP0, R1, k2Ref, dBPndC = " << _BP0 << _R1 << _k2Ref << _dBPndChallenge << _dBPndUptake;

    double toleranceCost = 0.01;
    double plasmaSearch = 1.;
    if ( level == 4 )
    {
        _parIncrements = {_BP0/5. , _R1/2. , _k2Ref/1., _tau4/1., _BP0/1. , plasmaSearch/1., _BP0/1.};
        toleranceCost = 0.05;
    }
    else if ( level == 3 )
    {
        _parIncrements = {_BP0/5. , _R1/4. , _k2Ref/2., _tau4/2., _BP0/2. , plasmaSearch/2., _BP0/2.};
        toleranceCost = 0.05;
    }
    else if ( level == 2 )
    {
        _parIncrements = {_BP0/10. , _R1/8. , _k2Ref/4., _tau4/4., _BP0/4. , plasmaSearch/4., _BP0/4.};
        toleranceCost = 0.05;
    }
    else if ( level == 1 )
    {
        _parIncrements = {_BP0/20. , _R1/16. , _k2Ref/8., _tau4/8., _BP0/8. , plasmaSearch/8., _BP0/8.};
        toleranceCost = 0.05;
    }

    int iCount = 0;
    // Find the best set of parameters at this resolution level.
    bool converged = false;
    int MAXIT = 50;

    dVector incPar = _parIncrements;
    while ( !converged && iCount < MAXIT )
    {
        double costInitial = runAndCalculateCost();
        ////////////////////////
        dVector costPar; costPar.fill(0.,_parValues.size());
        // Test each parameter by increasing and decreasing its value by the increment.
        for (int jPar=0; jPar<_parValues.size(); jPar++)
        {
            if ( _parAdjust[jPar] )
                // find the best value for this parameter.
                lineScan1D(jPar, costPar[jPar], incPar[jPar]);
            else
                costPar[jPar] = incPar[jPar]  = 0.;
        }

        // Normalize the costPar array.
        double mag=0.;
        for (int jPar=0; jPar<_parValues.size(); jPar++)
        {
            if ( _parAdjust[jPar] ) mag += SQR(costPar[jPar]);
        }
        mag = qSqrt(mag);
//        FUNC_INFO << "mag" << mag << "costPar" << costPar;
        // If the magnitude of the costPar vector is 0, then no values are lower than the
        // original value.
        if ( mag == 0. )
            break;
        for (int jPar=0; jPar<_parValues.size(); jPar++)
            costPar[jPar] /= mag;

//        FUNC_INFO << "incVector" << incPar;
//        FUNC_INFO << "costVector" << costPar;
        // Calculate the projection of the incPar along the normalized costPar vector.
        for (int jPar=0; jPar<_parValues.size(); jPar++)
        {
            incPar[jPar] = incPar[jPar] * costPar[jPar];
            if ( _parAdjust[jPar] ) _parValues[jPar] += incPar[jPar];
        }
//        FUNC_INFO << "_parValues result" << _parValues;
//        FUNC_INFO << "iteration" << iCount << "pars:" << _parValues;
//        FUNC_INFO << "iteration" << iCount << "incs:" << _parIncrements;
//        FUNC_INFO << "iteration" << iCount << "rCost:" << runAndCalculateCost()/costInitial;
        ////////////////////////
        double costFinal = runAndCalculateCost();
        converged = qAbs(costFinal/costInitial - 1.) < toleranceCost;
        converged &= iCount > 10;
//        FUNC_INFO << "converged?" << costInitial << costFinal << qAbs(costFinal/costInitial - 1.) << toleranceCost;
        iCount++;
    }
//    FUNC_EXIT << "BP0, R1, k2Ref = " << _BP0 << _R1 << _k2Ref << _dBPndChallenge << _dBPndUptake << "\n";
}

void simAnalyzer::lineScan1D( int iPar, double &costRelative, double &incrementOpt )
{ // return value is true if a change was made
//    FUNC_ENTER << iPar << costRelative << incrementOpt;
    costRelative = incrementOpt = 0.;
    if ( _parIncrements[iPar] <= 0 ) return;
//    FUNC_INFO << "value & increment1" << _parValues[iPar] << _parIncrements[iPar];
    double valueInitial = _parValues[iPar];
    double costInitial  = runAndCalculateCost();  // refit with all initial parameter set
    // Save the initial cost function.
    double xParabola[3], yParabola[3];
    xParabola[1] = valueInitial;
    yParabola[1] = costInitial;

    // Now test + and - directions
    // test the - direction
    double value0 = valueInitial - _parIncrements[iPar];
//    _parValues[iPar] = qMax(0.,value0);
    _parValues[iPar] = value0;
    double cost0 = runAndCalculateCost();
    xParabola[0] = value0;
    yParabola[0] = cost0;

    // test the + direction
    double value2 = valueInitial + _parIncrements[iPar];
//    _parValues[iPar] = qMax(0.,value2);
    _parValues[iPar] = value2;
    double cost2 = runAndCalculateCost();
    xParabola[2] = value2;
    yParabola[2] = cost2;
//    FUNC_INFO << "costs" << cost0 << cost2;

    // allocate values that might be required if we need to keep going in one direction to find the max
    bool capturedMinOrMax;
    if ( _optHigh )
        capturedMinOrMax = costInitial > cost0 && costInitial > cost2;  // costInitial is higher than either side, so good
    else
        capturedMinOrMax = costInitial < cost0 && costInitial < cost2;  // costInitial is lower than either side, so good
    if ( capturedMinOrMax )
    { // the increment range contains the maximum/minimum of the cost function, so interpolate to get the best estimate
        double xMax, yMax;
        if ( utilMath::ParabolicInterpolation(xParabola, yParabola, xMax, yMax) )
        { // set the final increment to the interpolated value; half the increment range for next time
            if ( _optHigh )
                costRelative = yMax - costInitial;  // cost function should be growing
            else
                costRelative = costInitial - yMax;  // cost function should be shrinking
//            FUNC_INFO << "xMax, yMax" << xMax << yMax;
            incrementOpt = xMax - valueInitial;
            _parIncrements[iPar] /= 2.;
        }
        else
            // This should never happen.
            costRelative = incrementOpt = 0.;
    }
    else if ( qFuzzyCompare(costInitial,cost0) || qFuzzyCompare(costInitial,cost2) )
        // not enough information
        costRelative = incrementOpt = 0.;
    else
    { // set the increment at the best edge
        if ( cost2 > cost0 )  // cost2 > costInitial > cost0
        {
            if ( _optHigh )
            {
                incrementOpt = _parIncrements[iPar];
                costRelative = cost2 - costInitial;
            }
            else
            {
                incrementOpt = - _parIncrements[iPar];
                costRelative = costInitial - cost0;
            }
        }
        else // cost0 > costInitial > cost2
        {
            if ( _optHigh )
            {
                incrementOpt = -_parIncrements[iPar];
                costRelative = cost0 - costInitial;
            }
            else
            {
                incrementOpt = _parIncrements[iPar];
                costRelative = costInitial - cost2;
            }
        }
    }
//    _parValues[iPar] = qMax(0.0,valueInitial);  // reset parameter value
    _parValues[iPar] = valueInitial;  // reset parameter value
//    FUNC_INFO << "iPar, value & increment2" << iPar << _parValues[iPar] << _parIncrements[iPar];
//    FUNC_EXIT << iPar << costRelative << incrementOpt;
    return;
}

void simAnalyzer::defineChallenge(int shape, double tau, dVector onTime, dVector offTime)
{
//    FUNC_ENTER;
    _challengeShape = shape;
    _challengeTau   = tau;
    _challengeOn    = onTime;
    _challengeOff   = offTime;
    updateChallengeShape();
}

void simAnalyzer::updateChallengeShape()
{
    _challengeCurve.fill(0.,_dtFine.size());
    if ( _challengeOn.size() == 0 ) return;

    for (int jt=0; jt<_dtFine.size(); jt++)
    {
        double time = getTimeFine(jt);
        for ( int jStim=0; jStim<_challengeOn.size(); jStim++ )
        {
            if ( _challengeShape == Challenge_Square )
            {
                if ( time >= _challengeOn[jStim] && time < _challengeOff[jStim] )
                    _challengeCurve[jt] += 1.;
            }
            else if ( _challengeShape == Challenge_RampUp )
            {
                double duration = _challengeOff[jStim] - _challengeOn[jStim];
                if ( time >= _challengeOn[jStim] && time < _challengeOff[jStim] )
                    _challengeCurve[jt] += ( time-_challengeOn[jStim] )/duration;
            }
            else if ( _challengeShape == Challenge_RampDown )
            {
                double duration = _challengeOff[jStim] - _challengeOn[jStim];
                if ( time >= _challengeOn[jStim] && time < _challengeOff[jStim] )
                    _challengeCurve[jt] += 1. - ( time-_challengeOn[jStim] )/duration;
            }
            else if ( _challengeShape == Challenge_Gamma )
            {
                double time0 = _challengeOn[jStim];
                double tau = _challengeTau;
                // double alpha = _challengeAlpha[indexChallenge];
                if ( time >= time0 )
                    _challengeCurve[jt] += (time-time0)/tau * qExp(1.-(time-time0)/tau);
            }
            else if ( _challengeShape == Challenge_Sigmoid )
            {
                double time0 = _challengeOn[jStim];
                double tau = _challengeTau;
                if ( time >= time0 )
                    _challengeCurve[jt] += (time-time0)/tau / qSqrt((1.+ (time-time0)/tau*(time-time0)/tau  ));
            }
        }
    }
    _challengeCurveBinned = downSample(_challengeCurve);
}
