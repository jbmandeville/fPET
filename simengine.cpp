#include "simengine.h"
#include <QtMath>
#include <QDebug>
#include <QRandomGenerator>

simEngine::simEngine()
{
    run();
}

void simEngine::run()
{
    generatePlasmaTAC();
    generateReferenceTAC();
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

    int nTime = static_cast<int>(_duration / _dt);
    double Cp_fast = 0.;  double Cp_slow = 0.;
    _Cp.clear();
    double magInfusion = 0.;
    if ( _KBol != 0. )
        magInfusion = _magBolus * exp(1) * _tauBolus / _KBol;
    for ( int jt=0; jt<nTime; jt++)
    {
        double time = jt * _dt;
        double plasmaInput = magInfusion + _magBolus * time/_tauBolus * qExp(1.-time/_tauBolus);
        if ( jt != 0. )
        {
            double dCpdt_fast  = _fracFast      * plasmaInput - _kFast * Cp_fast;
            double dCpdt_slow  = (1.-_fracFast) * plasmaInput - _kSlow * Cp_slow;
            Cp_fast += dCpdt_fast * _dt;
            Cp_slow += dCpdt_slow * _dt;
        }
        _Cp.append(Cp_fast + Cp_slow);
    }
    _CpBinned = downSample(_Cp);
}

void simEngine::generateReferenceTAC()
{
    // dCr_dt = K1 * Cp - k2 * Cr;
    int nTime = _Cp.size();
    double Cr = 0.;
    _Cr.clear();
    for ( int jt=0; jt<nTime; jt++)
    {
        double Cp = _Cp[jt];
        if ( jt != 0. )
        {
            double dCrdt = _K1Ref * Cp - _k2Ref * Cr;
            Cr += dCrdt * _dt;
        }
        _Cr.append(Cr + _fracPlasma * Cp);
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
    int nTime = _Cp.size();
    _Ct.clear();
    for ( int jt=0; jt<nTime; jt++)
    {
        int iDisplacementTime = static_cast<int>(_challengeTime / _dt);
        bool postChallenge = (jt > iDisplacementTime);
        if ( postChallenge && _deltaBPPercent != 0. )
        {
            double BP1 = _BP0 * (1. - _deltaBPPercent/100.);
            k2a = _k2 / (1. + BP1);
        }
        if ( jt > 0 )
        {
            double dCtdt = _K1 * _Cp[jt] - k2a * Ct ;
            Ct += dCtdt * _dt;
        }
        _Ct.append(Ct + _fracPlasma * _Cp[jt]);
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
    double Cf=0.;  double Cb=0.;
    int nTime = _Cp.size();
    _Cf.clear();  _Cb.clear();  _Ct.clear();
    for ( int jt=0; jt<nTime; jt++)
    {
        int iDisplacementTime = static_cast<int>(_challengeTime / _dt);
        bool postChallenge = (jt > iDisplacementTime);
        if ( postChallenge && _deltaBPPercent != 0. )
        {
            double BP1 = _BP0 * (1. - _deltaBPPercent/100.);
            k3 = BP1 * _k4;
        }
        if ( jt > 0 )
        {
            double dCfdt = _K1 * _Cp[jt] + _k4 * Cb - (_k2 + k3) * Cf;
            double dCbdt = k3 * Cf - _k4 * Cb;
            Cf += dCfdt * _dt;
            Cb += dCbdt * _dt;
        }
        _Cf.append(Cf);
        _Cb.append(Cb);
        _Ct.append(Cf + Cb + _fracPlasma * _Cp[jt]);
    }

    // Downsample the TAC and add Noise
    _CtBinned = downSample(_Ct);
    addNoise(_noiseTar, _CtBinned);
}

dVector simEngine::downSample(dVector original)
{
    dVector binned;
    // Downsample the TAC
    int nTime = original.size();
    for ( int jt=0; jt<nTime; jt+=_lDownSample)
    {
        double average = 0.;
        for (int jDown=0; jDown<_lDownSample; jDown++)
            average += original[jt+jDown];
        average /= static_cast<double>(_lDownSample);
        binned.append(average);
    }
    return binned;
}
void simEngine::addNoise(double noiseScale, dVector &timeVector)
{ // add noise to time vector AFTER down-sampling
    int nTime = timeVector.size();
    for ( int jt=0; jt<nTime; jt++)
    {
        // noise proportional to Cr * exp(t/T_1/2)
        double time = jt * _lDownSample * _dt;
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
      x = cutoff * ( -1. + 2. * QRandomGenerator::global()->generateDouble() );
      /* Compute the Gaussian function of x. */
      double arg = - 0.5 * (x*x) /sigma/sigma;
      yGauss = exp( arg );
      /* Choose a random y from 0 to 1. */
      y = QRandomGenerator::global()->generateDouble();
    }
  while (y > yGauss);

  return( x );
}
