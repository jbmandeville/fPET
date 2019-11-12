#include "simengine.h"
#include <QtMath>
#include <QDebug>
//#include <QRandomGenerator>

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
    qDebug() << "simEngine::initializeBins enter";
    if ( _nBins != _durationBin.size() )
        _durationBin.fill(1.,_nBins);  // duration in min
    if ( _nBins != _numberSamplesPerBin.size() )
        _numberSamplesPerBin.fill(10,_nBins);
    updateFineSamples();
}
void simEngine::updateFineSamples()
{
    _dtFine.clear();  _timeFine.clear(); _timeCoarse.clear();
    double time=0.;
    qDebug() << "simEngine::updateFineSamples bins" << _durationBin;
    for (int jBin=0; jBin<_nBins; jBin++)
    {
        double dtFine = _durationBin[jBin] / _numberSamplesPerBin[jBin];
        double avTime=0.;
        for (int jSample=0; jSample<_numberSamplesPerBin[jBin]; jSample++)
        {
            _dtFine.append(dtFine);
            _timeFine.append(time + dtFine/2.);
            avTime += _timeFine.last();
//            qDebug() << "avTime =" << avTime;
            time += dtFine;  // start of bin
        }
        avTime /= static_cast<double>(_numberSamplesPerBin[jBin]);
        _timeCoarse.append(avTime);
    }
//    qDebug() << "dtFine" << _dtFine;
//    qDebug() << "timeFine" << _timeFine;
//    qDebug() << "timeCoarse" << _timeCoarse;
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

    int nTime = _dtFine.size();
    double Cp_fast = 0.;  double Cp_slow = 0.;
    _Cp.clear();
    double magInfusion = 0.;
    if ( _KBol != 0. )
        magInfusion = _magBolus * exp(1) * _tauBolus / _KBol;
    for ( int jt=0; jt<nTime; jt++)
    {
        double time = _timeFine[jt];
        double plasmaInput = magInfusion + _magBolus * time/_tauBolus * qExp(1.-time/_tauBolus);
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
    // dCr_dt = K1 * Cp - k2 * Cr;
    int nTime = _dtFine.size();
//    qDebug() << "nTime =" << nTime << "dtFile =" << _dtFine;
    double Cr = 0.;
    _Cr.clear();
    for ( int jt=0; jt<nTime; jt++)
    {
        double Cp = _Cp[jt];
        if ( jt != 0. )
        {
            double dCrdt = _K1Ref * Cp - _k2Ref * Cr;
            Cr += dCrdt * _dtFine[jt];
        }
        // Cp is conc(tracer)/volume(blood), whereas Cr is conc(tracer)/volume(tissue), so multiple Cp by volume(blood)/volume(tissue)
        _Cr.append(Cr + _percentPlasmaRef/100. * Cp);
    }
//    qDebug() << "Cr =" << _Cr;

    // Downsample the TAC and add Noise
    _CrBinned = downSample(_Cr);
//    qDebug() << "CrBinnged =" << _CrBinned;
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
        if ( jt > 0 )
        {
            double dCtdt = _K1 * _Cp[jt] - k2a * Ct ;
            Ct += dCtdt * _dtFine[jt];
        }
        _Ct.append(Ct + _percentPlasmaTar/100. * _Cp[jt]);
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
        if ( jt > 0 )
        {
            double dCfdt = _K1 * _Cp[jt] + _k4 * Cb - (_k2 + k3) * Cf;
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
        yGauss = exp( arg );
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
        integral += tissue[jt] * _durationBin[jt];
    return integral;
}

dVector simEngine::differentiateTissueVector()
{
    int nTime = _CtBinned.size();
    dVector derivative;  derivative.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        double dt = _durationBin[jt];
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
        double dt = _durationBin[jt];
        double time = dt * jt;  // probably need to shift by 1/2 bin
        for ( int jtPrime=0; jtPrime<=jt; jtPrime++ )
        { // integrate from 0 to current time in run (=iTime_run & time)
            double timePrime = dt * jtPrime;
            convolution[jt] += tissue[jtPrime] * qExp(-(_k3+_k4)*(time-timePrime)) * dt;
        } // jtPrime
    } // jt
    return convolution;
}
double simEngine::getDuration()
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
    qDebug() << "compute dk2k3 truth" << _BP0 << BP1 << _k2 << k2k3_0 << k2k3_1 << dk2k3;
    return dk2k3/_k4;
    */
}
