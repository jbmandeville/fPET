#include "simengine.h"
#include <QtMath>
#include <QDebug>
//#include <QRandomGenerator> // include only in Qt5.10 and higher

#define RAND_0_1 (rand()/2147483647.)

simEngine::simEngine()
{
    // initialize
    initializeBins();
    // run
    run();
}
void simEngine::initializeBins()
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
void simEngine::updateFineSamples()
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

void simEngine::setCrFit(dVector CrFitCoarse)
{
    if ( _nBins != CrFitCoarse.size() )
    {
        FUNC_INFO << _nBins << CrFitCoarse.size();
        qFatal("Error in simEngine::setCrFit - the # of coarse bins should match the # in the fit vector");
    }
    _CrFitBinned = CrFitCoarse;

    tk::spline spline;
    spline.set_points(_timeCoarse,_CrFitBinned);    // currently it is required that X is already sorted

    // fine scale: interpolate LOESS using spline
    int nTimeFine = getNumberTimeBinsFine();
    _CrFit.clear();
    for (int jTime=0; jTime<nTimeFine; jTime++)
        _CrFit.append(qMax(0.,spline(_timeFine[jTime])));

    _CrFitDot.fill(0.,_CrFit.size());
    // Calculate the derivative
    for (int jBin=0; jBin<_CrFit.size(); jBin++)
    {
        if ( jBin != 0 )
            _CrFitDot[jBin] = (_CrFit[jBin] - _CrFit[jBin-1]) / _dtFine[jBin];
        else
            _CrFitDot[jBin] = 0.;
    }
}

void simEngine::run()
{
    if ( _simStartingPoint == simStart_fromPlasma )
    {
        generatePlasmaTAC();
        generateReferenceTAC();
    }
    generateTargetTAC();
}
void simEngine::generateTargetTAC()
{
    if ( _SRTM )
        generateTargetSRTM();
    else
        generateTargetFRTM();
}

void simEngine::generatePlasmaTAC()
{
    // For simplicity, treat the plasma as two separate compartments: fast and slow
    // dC_p/dt = I(t)-K_fast*C_p-_slow*C_p*+Ks_in*C_s*(V_p/V_s)

    int nTime = _dtFine.size();
    double Cp_fast = 0.;  double Cp_slow = 0.;
    _Cp.clear();
    double magInfusion = 0.;
    if ( _KBol != 0. )
        magInfusion = _magBolus * exp(1) * _tauBolus / _KBol;
    for ( int jt=0; jt<nTime; jt++)
    {
        double time = _timeFine[jt];
        double plasmaInput;
        if ( time < _KBolDelay )
            plasmaInput = 0.          + _magBolus * time/_tauBolus * qExp(1.-time/_tauBolus);
        else
            plasmaInput = magInfusion + _magBolus * time/_tauBolus * qExp(1.-time/_tauBolus);
        if ( jt != 0. )
        {
            double dCpdt_fast  = _fracFast      * plasmaInput - _kFast * Cp_fast;
            double dCpdt_slow  = (1.-_fracFast) * plasmaInput - _kSlow * Cp_slow;
            Cp_fast += dCpdt_fast * _dtFine[jt];
            Cp_slow += dCpdt_slow * _dtFine[jt];
        }
        _Cp.append(Cp_fast + Cp_slow);
    }
    _CpBinned = downSample(_Cp);
}

double simEngine::baselineBasisFunction(int iPoly, double x)
{ // return Legendre polynomial of x with order iPoly, where x is symmetric in (-1.,1.); x = -1. + 2.*time/duration
    double value;
    if ( iPoly == 0 )
        value = 1.;
    else if ( iPoly == 1 )
        value = x;
    else if ( iPoly == 2 )
        value = 0.5 * (3*x*x -1.);
    else if ( iPoly == 3 )
        value = 0.5 * (5.*x*x*x - 3.*x);
    else if ( iPoly == 4 )
        value = 0.125 * (35.*x*x*x*x - 30.*x*x + 3.);
    else // if ( iPoly == 5 )
        value = 0.125 * (63.*x*x*x*x*x - 70.*x*x*x + 15.*x);
    return value;
}
double simEngine::gammaVariateFunction(double time, double onset, double alpha, double tau)
{ // return a normalized gamma variate function of form t^alpha * exp(alpha*t)
    double t = (time - onset) / tau;
    return qPow(t,alpha) * qExp(alpha*(1.-t));
}

void simEngine::generateReferenceTAC()
{
    // dCr_dt = K1 * Cp - k2 * Cr;
    // Cp = ( dCr_dt + k2 * Cr) / K1;

    int nTime = _dtFine.size();
    double Cr = 0.;
    _Cr.clear();  _CrDot.clear();
    for ( int jt=0; jt<nTime; jt++)
    {
        double Cp = _Cp[jt];
        if ( jt != 0. )
        {
            double dCrdt = _K1Ref * Cp - _k2Ref * Cr;
            Cr += dCrdt * _dtFine[jt];
            _CrDot.append(dCrdt);
        }
        // Cp is conc(tracer)/volume(blood), whereas Cr is conc(tracer)/volume(tissue), so multiple Cp by volume(blood)/volume(tissue)
        _Cr.append(Cr + _percentPlasmaRef/100. * Cp);
    }

    // Downsample the TAC and add Noise
    _CrBinned = downSample(_Cr);
    addNoise(_noiseRef, _CrBinned);
}

void simEngine::generateTargetSRTM()
{
    // dCt_dt = K1 * Cp - k2a * Ct

    // Derived quantities: update tissue properties
    _K1  = _R1 * _K1Ref;
    _k2  = _R1 * _k2Ref;
    double k2a = _k2 / (1. + _BP0);

    double Ct=0.;
    int nTime = _dtFine.size();
    _Ct.clear();
    for ( int jt=0; jt<nTime; jt++)
    {
        int iDisplacementTime = static_cast<int>(_challengeTime / _dtFine[jt]);
        bool postChallenge = (jt > iDisplacementTime);
        if ( postChallenge && _deltaBPPercent != 0. )
        {
            double BP1 = _BP0 * (1. - _deltaBPPercent/100.);
            k2a = _k2 / (1. + BP1);
        }
        dVector Cp;
        if ( _simStartingPoint == simStart_fromPlasma )
        {
            FUNC_INFO << "start from plasma";
            Cp = _Cp;
        }
        else
        {
            FUNC_INFO << "start from fit";
            int nTime = _CrFit.size();
            Cp.resize(nTime);
            for ( int jTime=0; jTime<nTime; jTime++ )
                Cp[jTime] = _k2Ref/_K1Ref * _CrFit[jTime] + _CrFitDot[jTime]/_K1Ref;
        }
        if ( jt > 0 )
        {
            double dCtdt = _K1 * Cp[jt] - k2a * Ct ;
            Ct += dCtdt * _dtFine[jt];
        }
        _Ct.append(Ct + _percentPlasmaTar/100. * Cp[jt]);
    }

    // Downsample the TAC and add Noise
    _CtBinned = downSample(_Ct);
    addNoise(_noiseTar, _CtBinned);
}

void simEngine::generateTargetFRTM()
{
    // dCb_dt = k3 * Cf - koff * Cb = kon * Bavail * Cf - koff * Cb
    // dCf_dt = K1 * Cp + koff * Cb - (k2+k3)*Cf

    // Derived quantities: update tissue properties
    _K1  = _R1 * _K1Ref;
    _k2  = _R1 * _k2Ref;
    _k3  = _BP0  * _k4;

    double k3=_k3;
    double Cf=0.;  double Cb=0.;  double lastCr=0.;
    int nTime = _dtFine.size();
    _Cf.clear();  _Cb.clear();  _Ct.clear();
    for ( int jt=0; jt<nTime; jt++)
    {
        int iDisplacementTime = static_cast<int>(_challengeTime / _dtFine[jt]);
        bool postChallenge = (jt > iDisplacementTime);
        if ( postChallenge && _deltaBPPercent != 0. )
        {
            double BP1 = _BP0 * (1. - _deltaBPPercent/100.);
            k3 = BP1 * _k4;
        }
        double Cr = _CrFit[jt];
        if ( jt > 0 )
        {
            double dCrdt = (Cr - lastCr)/_dtFine[jt];
            double dCfdt;
            if ( _simStartingPoint == simStart_fromPlasma )
                dCfdt = _K1 * _Cp[jt] + _k4 * Cb - (_k2 + k3) * Cf;
            else
                dCfdt = _K1/_K1Ref * ( _k2Ref * _CrFit[jt] + dCrdt) + _k4 * Cb - (_k2 + k3) * Cf;
            double dCbdt = k3 * Cf - _k4 * Cb;
            Cf += dCfdt * _dtFine[jt];
            Cb += dCbdt * _dtFine[jt];
            lastCr = Cr;
        }
        _Cf.append(Cf);
        _Cb.append(Cb);
        _Ct.append(Cf + Cb + _percentPlasmaTar/100. * _Cp[jt]);
        lastCr = Cr;
    }

    // Downsample the TAC and add Noise
    _CtBinned = downSample(_Ct);
    addNoise(_noiseTar, _CtBinned);
}

dVector simEngine::downSample(dVector original)
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
void simEngine::addNoise(double noiseScale, dVector &timeVector)
{ // add noise to time vector AFTER down-sampling
    int nTime = timeVector.size();
    for ( int jt=0; jt<nTime; jt++)
    {
        // noise proportional to Cr * exp(t/T_1/2)
        double time = _timeCoarse[jt];
        double noiseMag = qSqrt(noiseScale * timeVector[jt] * exp(time/28.9));
        double noise = GaussianRandomizer(noiseMag, 2.*noiseMag);
        timeVector[jt] += noise;
    }
}
double simEngine::GaussianRandomizer(double sigma, double cutoff)
{
    double yGauss, y, x;
    do
    {
        /* Choose x between +- the cutoff */
//        x = cutoff * ( -1. + 2. * QRandomGenerator::global()->generateDouble() );
        x = cutoff * ( -1. + 2. * RAND_0_1 );
        /* Compute the Gaussian function of x. */
        double arg = - 0.5 * (x*x) /sigma/sigma;
        yGauss = qExp( arg );
        /* Choose a random y from 0 to 1. */
//        y = QRandomGenerator::global()->generateDouble();
        y = RAND_0_1;
    }
    while (y > yGauss);

  return( x );
}

dVector simEngine::FRTMNewAnalyticalSolution()
{
    int nTime = _CpBinned.size();
    dVector convolution = calculateConvolution(_CtBinned);

    dVector Ct; Ct.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
        Ct[jt] = _R1 * _CrBinned[jt] - _k2 * (integralOf(_CtBinned,jt) - integralOf(_CrBinned,jt)) + _k2*_k3 * integralOf(convolution,jt);
    return Ct;
}

dVector simEngine::FRTMOldAnalyticalSolution()
{
    int nTime = _CpBinned.size();

    // Calculate the derivative of the tissue vector
    dVector tissDerivative = differentiateTissueVector();
    dVector convolution = calculateConvolution(tissDerivative);
    double k2a = _k2 / (1. + _k3/_k4);

    dVector CrTerm; CrTerm.resize(nTime);
    dVector CtTerm; CtTerm.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        CrTerm[jt] = _CrBinned[jt]-convolution[jt];
        CtTerm[jt] = _CtBinned[jt]-convolution[jt];
    }

    dVector Ct; Ct.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
        Ct[jt] = _R1 * _CrBinned[jt] + _k2 * integralOf(CrTerm,jt) - k2a * integralOf(CtTerm,jt);
    return Ct;
}

double simEngine::integralOf(dVector tissue, int iTime)
{
    double integral=0;
    for ( int jt=0; jt<=iTime; jt++)
    {
        double dt = static_cast<double>(_durationBinSec[jt])/60.;
        integral += tissue[jt] * dt;
    }
    return integral;
}

dVector simEngine::differentiateTissueVector()
{
    int nTime = _CtBinned.size();
    dVector derivative;  derivative.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        double dt = static_cast<double>(_durationBinSec[jt])/60.;
        if ( jt != 0. )
            derivative[jt] = (_CtBinned[jt] - _CtBinned[jt-1]) / dt;
        else
            derivative[jt] = _CtBinned[jt] / dt;
    }
    return derivative;
}
dVector simEngine::calculateConvolution(dVector tissue)
{
    int nTime = tissue.size();
    dVector convolution;  convolution.fill(0.,nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        double dt = static_cast<double>(_durationBinSec[jt])/60.;
        double time = dt * jt;  // probably need to shift by 1/2 bin
        for ( int jtPrime=0; jtPrime<=jt; jtPrime++ )
        { // integrate from 0 to current time in run (=iTime_run & time)
            double timePrime = dt * jtPrime;
            convolution[jt] += tissue[jtPrime] * qExp(-(_k3+_k4)*(time-timePrime)) * dt;
        } // jtPrime
    } // jt
    return convolution;
}
double simEngine::getDurationScan()
{
    double duration=0.;
    for (int jt=0; jt<_dtFine.size(); jt++)
        duration += _dtFine[jt];
    return duration;
}
double simEngine::getdk2a()
{
    // k2a = k2 / (1+_BP0)
    // BP1 = BP0
    double k2a_0 = _k2/(1.+_BP0);
    double BP1 = _BP0 * (1. - _deltaBPPercent/100.);
    double k2a_1 = _k2/(1.+BP1);
    double dk2a = k2a_1 - k2a_0;
    return dk2a;
}
double simEngine::getdk2k3()
{
    double BP1 = _BP0 * (1. - _deltaBPPercent/100.);
    double DBP = _BP0 - BP1;
    return - _k2 * DBP;  // don't use PET convention
    /*
    double k2k3_0 = _k2 * _BP0;
    double k2k3_1 = _k2 * BP1;
    double dk2k3 = k2k3_0 - k2k3_1;
    return dk2k3/_k4;
    */
}
