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
    _parValuesFRTM3 = {_BP0,_R1,_k2Ref};
    runAndCalculateCost();
}
void simAnalyzer::initializeBins()
{
    if ( _nBins != _durationBinSec.size() )
    {
        double defaultBinSize=60;
        if ( _durationBinSec.size() > 0 ) defaultBinSize = _durationBinSec[0];
        _durationBinSec.fill(60,_nBins);  // duration in min
    }
    if ( _nBins != _numberSamplesPerBin.size() )
        _numberSamplesPerBin.fill(10,_nBins);
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
    FUNC_ENTER << _parValuesFRTM3;
    _BP0   = _parValuesFRTM3.at(0);
    _R1    = _parValuesFRTM3.at(1);
    _k2Ref = _parValuesFRTM3.at(2);
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
        //
        _sigma2 += (_CtData[jt]-_CtBinned[jt]) * (_CtData[jt]-_CtBinned[jt]);
//        }
    }
    double dof = static_cast<double>(nTime);
    _sigma2 /= dof;
    return _sigma2;
}

void simAnalyzer::generateTargetTACFRTM()
{
    // dCr_dt = K1Ref * Cp - k2Ref * Cr
    // dCb_dt = k3 * Cf - koff * Cb = kon * Bavail * Cf - koff * Cb
    // dCf_dt = K1 * Cp + koff * Cb - (k2+k3)*Cf
    // R1 = K1/K1Ref = k2/k2Ref

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
        int iDisplacementTime = static_cast<int>(_challengeTime / _dtFine[jt]);
        bool postChallenge = (jt > iDisplacementTime);
        if ( postChallenge && _deltaBPPercent != 0. )
        {
            double BP1 = _BP0 * (1. - _deltaBPPercent/100.);
            k3 = BP1 * k4;
        }
        double Cr = _Cr[jt];
        // To add plasma contamination, one needs K1Rref (not just R1)
        if ( jt > 0 )
        {
            double dCrdt = (Cr - lastCr)/_dtFine[jt];
            double dCfdt;
//            Cp = ( _k2Ref * _Cr[jt] + dCrdt) / _K1Ref;        // plasma
//            dCfdt = _K1 * Cp + k4 * Cb - (_k2 + k3) * Cf;    // dCfdt from plasma

//            dCfdt = _R1 * ( _k2Ref * _Cr[jt] + dCrdt) + k4 * Cb - (_k2 + k3) * Cf;
            dCfdt = _R1 * ( _k2Ref * _Cr[jt] + _CrDot[jt]) + k4 * Cb - (_k2 + k3) * Cf;
            double dCbdt = k3 * Cf - k4 * Cb;
            Cf += dCfdt * _dtFine[jt];
            Cb += dCbdt * _dtFine[jt];
            lastCr = Cr;
        }
//        _Ct.append(Cf + Cb + _percentPlasmaTar/100. * Cp);
        _Ct.append(Cf + Cb);
        lastCr = Cr;
    }
    // Downsample the TAC and add Noise
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
        int iDisplacementTime = static_cast<int>(_challengeTime / _dtFine[jt]);
        bool postChallenge = (jt > iDisplacementTime);
        if ( postChallenge && _deltaBPPercent != 0. )
        {
            double BP1 = _BP0 * (1. - _deltaBPPercent/100.);
            k2a = _k2 / (1. + BP1);
        }
        double Cr = _Cr[jt];
        // To add plasma contamination, one needs K1Rref (not just R1)
        if ( jt > 0 )
        {
            double dCrdt = (Cr - lastCr)/_dtFine[jt];
//            Cp = ( _k2Ref * _Cr[jt] + dCrdt) / _K1Ref;        // plasma
            double dCtdt = _R1 * ( _k2Ref * _Cr[jt] + dCrdt) - k2a * Ct;
            Ct += dCtdt * _dtFine[jt];
            lastCr = Cr;
        }
//        _Ct.append(Cf + Cb + _percentPlasmaTar/100. * Cp);
        _Ct.append(Ct);
        lastCr = Cr;
    }
    // Downsample the TAC and add Noise
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
    double BP1 = _BP0 * (1. - _deltaBPPercent/100.);
    double k2a_1 = _k2/(1.+BP1);
    double dk2a = k2a_1 - k2a_0;
    return dk2a;
}
double simAnalyzer::getdk2k3()
{
    double BP1 = _BP0 * (1. - _deltaBPPercent/100.);
    double DBP = _BP0 - BP1;
    return - _k2 * DBP;  // don't use PET convention
}

void simAnalyzer::fit(double BP0, double R1, double k2Ref, bool fitAll)
{
    FUNC_ENTER << "BP0, R1, k2Ref = " << BP0 << R1 << k2Ref;
    _BP0   = BP0;
    _R1    = R1;
    _k2Ref = k2Ref;
    if ( fitAll )
        _parAdjust = {true, true, true};
    else
        _parAdjust = {true, true, false};
    fitTACByGridSearch(6);  // 6 -> 1
    fitTACByLineScan(4);
    fitTACByLineScan(4);
    fitTACByLineScan(3);
    fitTACByLineScan(3);
    fitTACByLineScan(2);
    fitTACByLineScan(2);
    fitTACByLineScan(1);
}

void simAnalyzer::fitTACByGridSearch(int level)
{
//    FUNC_ENTER << "BP0, R1, k2Ref = " << _BP0 << _R1 << _k2Ref << "level" << level;
    // Fit a gamma function plus a polynomial by repeated GLM using a grid search for the best gamma parameters (3)
    double BP0   = _BP0;
    double R1    = _R1;
    double k2Ref = _k2Ref;

    double BP0Width, R1Width, k2RefWidth;
    if ( level == 6 )
    {
        BP0Width = BP0/5.;
        k2RefWidth = k2Ref/1.;
        R1Width  = R1/2.;
    }
    else if ( level == 5 )
    {
        BP0Width = BP0/10.;
        k2RefWidth = k2Ref/2.;
        R1Width  = R1/4.;
    }
    else if ( level == 4 )
    {
        BP0Width = BP0/20.;
        k2RefWidth = k2Ref/4.;
        R1Width  = R1/8.;
    }
    else if ( level == 3 )
    {
        BP0Width = BP0/40.;
        k2RefWidth = k2Ref/8.;
        R1Width  = R1/16.;
    }
    else if ( level == 2 )
    {
        BP0Width = BP0/80.;
        k2RefWidth = k2Ref/16.;
        R1Width  = R1/32.;
    }
    else // if ( level == 1 )
    {
        BP0Width = BP0 * 160.;
        k2RefWidth = k2Ref/32.;
        R1Width  = R1/64.;
    }

//    FUNC_INFO << "widths for BP0Width, R1Width, k2RefWidth" << BP0Width << R1Width << k2RefWidth;

    double stepDivider = 5.;
    double BP0Step   = BP0Width  /stepDivider;
    double R1Step    = R1Width   /stepDivider;
    double k2RefStep = k2RefWidth/stepDivider;

    double BP0Start   = qMax(0.,BP0-BP0Width);      double BP0Stop = BP0 + BP0Width;
    double R1Start    = qMax(0.,R1-R1Width);        double R1Stop   = R1 + R1Width;
    double k2RefStart = qMax(0.,k2Ref-k2RefWidth);  double k2RefStop = k2Ref + k2RefWidth;

    _parValuesFRTM3     = {_BP0,_R1,_k2Ref};
    _parIncrementsFRTM3 = {BP0Step, R1Step, k2RefStep};
    dVector bestSet = _parValuesFRTM3;
    double bestSigma2 = runAndCalculateCost();

    if ( !_parAdjust[0] )
        BP0Start = BP0Step = _BP0;
    if ( !_parAdjust[1] )
        R1Start = R1Stop = _R1;
    if ( !_parAdjust[2] )
        k2RefStart = k2RefStop = _k2Ref;

    for ( BP0=BP0Start; BP0<BP0Stop; BP0+=BP0Step)
    {
        _BP0 = BP0;
        for ( R1=R1Start; R1<R1Stop; R1+=R1Step)
        {
            _R1 = R1;
            for ( k2Ref=k2RefStart; k2Ref<k2RefStop; k2Ref+=k2RefStep)
            {
                _k2Ref = k2Ref;
                _parValuesFRTM3 = {_BP0,_R1,_k2Ref};  double sigma2 = runAndCalculateCost();
                if ( sigma2 < bestSigma2 )
                {
                    bestSigma2 = sigma2;
                    bestSet = _parValuesFRTM3;
                }
            }
        }
    }
    _parValuesFRTM3 = bestSet;
    bestSigma2 = runAndCalculateCost();
//    FUNC_INFO << "output " << bestSet << "with sigma =" << qSqrt(bestSigma2) << "and tau4 = " << _tau4 << "\n";
}

void simAnalyzer::fitTACByLineScan(int level) // toleranceCost is a fraction (e.g., 0.01 is 1%)
{
    FUNC_ENTER;

    double toleranceCost = 0.01;
    if ( level == 4 )
    {
        _parIncrementsFRTM3 = {_BP0/5. , _R1/2. , _k2Ref/1.};
        toleranceCost = 0.05;
    }
    else if ( level == 3 )
    {
        _parIncrementsFRTM3 = {_BP0/5. , _R1/4. , _k2Ref/2.};
        toleranceCost = 0.05;
    }
    else if ( level == 2 )
    {
        _parIncrementsFRTM3 = {_BP0/10. , _R1/8. , _k2Ref/4.};
        toleranceCost = 0.05;
    }
    else if ( level == 1 )
    {
        _parIncrementsFRTM3 = {_BP0/20. , _R1/16. , _k2Ref/8.};
        toleranceCost = 0.05;
    }

    int iCount = 0;
    // Find the best set of parameters at this resolution level.
    bool converged = false;
    int MAXIT = 50;

    dVector incPar = _parIncrementsFRTM3;
    while ( !converged && iCount < MAXIT )
    {
        double costInitial = runAndCalculateCost();
        ////////////////////////
        dVector costPar; costPar.fill(0.,_parValuesFRTM3.size());
        // Test each parameter by increasing and decreasing its value by the increment.
        for (int jPar=0; jPar<_parValuesFRTM3.size(); jPar++)
        {
            if ( _parAdjust[jPar] )
                // find the best value for this parameter.
                lineScan1D(jPar, costPar[jPar], incPar[jPar]);
            else
                costPar[jPar] = incPar[jPar]  = 0.;
        }

        // Normalize the costPar array.
        double mag=0.;
        for (int jPar=0; jPar<_parValuesFRTM3.size(); jPar++)
            mag += SQR(costPar[jPar]);
        mag = qSqrt(mag);
        FUNC_INFO << "mag" << mag << "costPar" << costPar;
        // If the magnitude of the costPar vector is 0, then no values are lower than the
        // original value.
        if ( mag == 0. )
            break;
        for (int jPar=0; jPar<_parValuesFRTM3.size(); jPar++)
            costPar[jPar] /= mag;

        FUNC_INFO << "incVector" << incPar;
        FUNC_INFO << "costVector" << costPar;
        // Calculate the projection of the incPar along the normalized costPar vector.
        for (int jPar=0; jPar<_parValuesFRTM3.size(); jPar++)
        {
            incPar[jPar] = incPar[jPar] * costPar[jPar];
            _parValuesFRTM3[jPar] += incPar[jPar];
            _parValuesFRTM3[jPar] = qMax(0.,_parValuesFRTM3[jPar]);
        }
        FUNC_INFO << "_parValuesFRTM3 result" << _parValuesFRTM3;
        FUNC_INFO << "iteration" << iCount << "pars:" << _parValuesFRTM3[0] << _parValuesFRTM3[1] << _parValuesFRTM3[2];
        FUNC_INFO << "iteration" << iCount << "incs:" << _parIncrementsFRTM3[0] << _parIncrementsFRTM3[1] << _parIncrementsFRTM3[2];
        FUNC_INFO << "iteration" << iCount << "rCost:" << runAndCalculateCost()/costInitial;
        ////////////////////////
        double costFinal = runAndCalculateCost();
        converged = qAbs(costFinal/costInitial - 1.) < toleranceCost;
        converged &= iCount > 10;
        FUNC_INFO << "converged?" << costInitial << costFinal << qAbs(costFinal/costInitial - 1.) << toleranceCost;
        iCount++;
    }
    FUNC_ENTER << "\n";
}

void simAnalyzer::lineScan1D( int iPar, double &costRelative, double &incrementOpt )
{ // return value is true if a change was made
//    FUNC_ENTER << iPar << costRelative << incrementOpt;
//    FUNC_INFO << "value & increment1" << _parValuesFRTM3[iPar] << _parIncrementsFRTM3[iPar];
    double valueInitial = _parValuesFRTM3[iPar];
    double costInitial  = runAndCalculateCost();  // refit with all initial parameter set
    // Save the initial cost function.
    double xParabola[3], yParabola[3];
    xParabola[1] = valueInitial;
    yParabola[1] = costInitial;

    // Now test + and - directions
    // test the - direction
    double value0 = valueInitial - _parIncrementsFRTM3[iPar];
    _parValuesFRTM3[iPar] = qMax(0.01,value0);
    double cost0 = runAndCalculateCost();
    xParabola[0] = value0;
    yParabola[0] = cost0;

    // test the + direction
    double value2 = valueInitial + _parIncrementsFRTM3[iPar];
    _parValuesFRTM3[iPar] = qMax(0.01,value2);
    double cost2 = runAndCalculateCost();
    xParabola[2] = value2;
    yParabola[2] = cost2;

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
            incrementOpt = xMax - valueInitial;
            _parIncrementsFRTM3[iPar] /= 2.;
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
                incrementOpt = _parIncrementsFRTM3[iPar];
                costRelative = cost2 - costInitial;
            }
            else
            {
                incrementOpt = - _parIncrementsFRTM3[iPar];
                costRelative = costInitial - cost0;
            }
        }
        else // cost0 > costInitial > cost2
        {
            if ( _optHigh )
            {
                incrementOpt = -_parIncrementsFRTM3[iPar];
                costRelative = cost0 - costInitial;
            }
            else
            {
                incrementOpt = _parIncrementsFRTM3[iPar];
                costRelative = costInitial - cost2;
            }
        }
    }
    _parValuesFRTM3[iPar] = qMax(0.01,valueInitial);  // reset parameter value
//    FUNC_INFO << "value & increment2" << _parValuesFRTM3[iPar] << _parIncrementsFRTM3[iPar];
//    FUNC_EXIT << iPar << costRelative << incrementOpt;
    return;
}

