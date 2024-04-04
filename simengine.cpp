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

void simEngine::run()
{
    FUNC_ENTER;
    generatePlasmaTAC();
    generateTargetTAC();
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
        double time = _timeFine[jt] - _timeCoarse[0];
        double plasmaInput = 0.;
        if ( time > 0. )
        {
            plasmaInput = _magBolus * time/_tauBolus * qExp(1.-time/_tauBolus);
            if ( time > _KBolDelay ) plasmaInput += magInfusion;
            double dCpdt_fast  = _fracFast      * plasmaInput - _kFast * Cp_fast;
            double dCpdt_slow  = (1.-_fracFast) * plasmaInput - _kSlow * Cp_slow;
            Cp_fast += dCpdt_fast * _dtFine[jt];
            Cp_slow += dCpdt_slow * _dtFine[jt];
            _Cp.append(Cp_fast + Cp_slow);
        }
        else
            _Cp.append(0.);
    }
    _CpBinned = downSample(_Cp);
}

void simEngine::generateTargetTAC()
{
    FUNC_ENTER;
    // dCf_dt(t) = K1 * Cp(t) - [k2 + k3(t)] * Cf(t) = K1 * Cp(t) - k3(0)*[1+k2/k3(0)+dk3(t)/k3(0)] * Cf(t)
    // dCb_dt = k3(t) * Cf = k3(0) * [1 + dk3(t)/k3(0)] * Cf
    // parameters:
    // K1, k3(0), k2/k3(0), dk3(t)/k3(0)

    // Derived quantities: update tissue properties
    double k2  = _k2_div_k3 * _k3;

    double k3;
    double Cf=0.;  double Cb=0.;
    int nTime = _dtFine.size();
    _Cf.clear();  _Cb.clear();  _Ct.clear();
    for ( int jt=0; jt<nTime; jt++)
    {
        FUNC_INFO << "jt" << jt;
        int iOnset1  = static_cast<int>(_challengeOnsetTime1  / _dtFine[jt]);
        int iOffset1 = static_cast<int>(_challengeOffsetTime1 / _dtFine[jt]);
        int iOnset2  = static_cast<int>(_challengeOnsetTime2  / _dtFine[jt]);
        int iOffset2 = static_cast<int>(_challengeOffsetTime2 / _dtFine[jt]);
        bool duringChallenge1 = (jt > iOnset1 && jt < iOffset1);
        bool duringChallenge2 = (jt > iOnset2 && jt < iOffset2);
        if ( duringChallenge1 && _dk3_percent1 != 0. )
        {
            double k3total = _k3 * (1. + _dk3_percent1/100.);
            k3 = k3total;
        }
        else if ( duringChallenge2 && _dk3_percent2 != 0. )
        {
            double k3total = _k3 * (1. + _dk3_percent2/100.);
            k3 = k3total;
        }
        else
            k3 = _k3;
        if ( jt > 0 )
        {
            FUNC_INFO << "k3[" << jt << "] =" << k3;
            double dCfdt = _K1 * _Cp[jt] - (k2 + k3) * Cf;
            double dCbdt = k3 * Cf;
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
//    double decayTimeConstant = 28.9;
    double decayTimeConstant = 151.;
    int nTime = timeVector.size();
    for ( int jt=0; jt<nTime; jt++)
    {
        // noise proportional to Cr * exp(t/T_1/2)
        double time = _timeCoarse[jt];
        double noiseMag = qSqrt(noiseScale * timeVector[jt] * exp(time/decayTimeConstant));
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
