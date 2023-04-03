#include "simengine.h"
#include <QtMath>
#include <QDebug>
//#include <QRandomGenerator> // include only in Qt5.10 and higher

#define RAND_0_1 (rand()/2147483647.)

simEngine::simEngine()
{
    FUNC_ENTER;
    // initialize
    initializeBins();
    // run
    run();
}
void simEngine::initializeBins()
{
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
    if ( _simStartingPoint == simStart_fromPlasma )
    {
        generatePlasmaTAC();
        generateReferenceTAC();
    }
    generateTargetTAC();
}
void simEngine::generateTargetTAC()
{
    FUNC_ENTER;
    if ( _SRTM )
        generateTargetSRTM();
    else
        generateTargetFRTM();
}

void simEngine::generatePlasmaTAC()
{
    FUNC_ENTER;
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

void simEngine::generateReferenceTAC()
{
    FUNC_ENTER;
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
    FUNC_ENTER;
    // dCt_dt = K1 * Cp - k2a * Ct

    // Derived quantities: update tissue properties
    _K1  = _R1 * _K1Ref;
    double k2  = _R1 * _k2Ref / _DV;
    double k2a = k2 / (1. + _BP0);

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
            k2a = k2 / (1. + BP1);
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
    FUNC_ENTER;
    // dCb_dt = k3 * Cf - koff * Cb = kon * Bavail * Cf - koff * Cb
    // dCf_dt = K1 * Cp + koff * Cb - (k2+k3)*Cf

    // Derived quantities: update tissue properties
    _K1  = _R1 * _K1Ref;
    double k2  = _R1 * _k2Ref / _DV;
    _k3  = _BP0  * _k4;

    double k3=_k3;
    double Cf=0.;  double Cb=0.;  double lastCr=0.;
    int nTime = _dtFine.size();
    _Cf.clear();  _Cb.clear();  _Ct.clear();
    FUNC_INFO << "sizes" << _dtFine.size() << _CrFit.size();
    FUNC_INFO << "_simStartingPoint" << _simStartingPoint;
    for ( int jt=0; jt<nTime; jt++)
    {
        FUNC_INFO << "jt" << jt;
        int iDisplacementTime = static_cast<int>(_challengeTime / _dtFine[jt]);
        bool postChallenge = (jt > iDisplacementTime);
        if ( postChallenge && _deltaBPPercent != 0. )
        {
            double BP1 = _BP0 * (1. - _deltaBPPercent/100.);
            k3 = BP1 * _k4;
        }
        if ( jt > 0 )
        {
            double dCfdt;
            if ( _simStartingPoint == simStart_fromPlasma )
                dCfdt = _K1 * _Cp[jt] + _k4 * Cb - (k2 + k3) * Cf;
            else
            {
                double Cr = _CrFit[jt];
                double dCrdt = (Cr - lastCr)/_dtFine[jt];
                dCfdt = _K1/_K1Ref * ( _k2Ref * Cr + dCrdt) + _k4 * Cb - (k2 + k3) * Cf;
                lastCr = Cr;
            }
            double dCbdt = k3 * Cf - _k4 * Cb;
            Cf += dCfdt * _dtFine[jt];
            Cb += dCbdt * _dtFine[jt];
        }
        _Cf.append(Cf);
        _Cb.append(Cb);
        _Ct.append(Cf + Cb + _percentPlasmaTar/100. * _Cp[jt]);
    }

    // Downsample the TAC and add Noise
    _CtBinned = downSample(_Ct);
    addNoise(_noiseTar, _CtBinned);
    FUNC_EXIT;
}

dVector simEngine::downSample(dVector original)
{
    FUNC_ENTER;
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
    FUNC_ENTER;
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
    FUNC_ENTER;
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
