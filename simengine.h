#ifndef SIMENGINE_H
#define SIMENGINE_H

#include <QDebug>

#include "io.h"
#include "generalglm.h"

class simEngine
{
private:
    // Input quantities ///////////////////////
    // Setup bins
    int _nBins = 90;
    iVector _durationBinSec;     // [_nBins]; seconds
    iVector _numberSamplesPerBin;  // [_nBins]
    dVector _dtFine;          // [sum of _numberSamplesPerBin = # of fine bins]
    dVector _timeFine;        // [sum of _numberSamplesPerBin = # of fine bins]
    dVector _dtCoarse;        // [_nBins]
    dVector _timeCoarse;      // [_nBins]
    // Infusion
    double _magBolus=15.e5;
    double _tauBolus=0.25;
    double _KBol=45.;
    // Elimination
    double _fracFast=0.9;   // fraction of elimination that is fast
    double _kFast=1./1.2;    // 1/T_alpha fast elimination
    double _kSlow=1./20.;    // 1/T_beta  slow elimination
    double _percentPlasmaRef=0.0 ; // CBV fraction (set to 0 to ignore plasma)
    double _percentPlasmaTar=0.0 ; // CBV fraction (set to 0 to ignore plasma)
    // Reference region
    double _K1Ref=1./15.;
    double _k2Ref=1./2.5;
    // Target region: state
    double _k4=1./10.;
    double _R1=1.;
    double _BP0=4.;     // BP = k3/koff = kon*Bavail0/koff
    // Challenge
    double _challengeTime=40.;
    double _deltaBPPercent=0.;  // deltaBPnd in %
    bool _SRTM = false;  // true if 1/k4 is set to zero (k4 = infinity)

    // Where to start
    bool _startWithPlasma=true;
    // structure for fitting the RR
    GeneralGLM _glmRR;

    // derived quantities //////////////////////////
    double _K1;    // K1 = R1 * K1_ref
    double _k2;    // k2 = R1 * k2_ref
    double _k3;    // k3 = BP * koff
    // noise
    double _noiseRef=0.;
    double _noiseTar=0.;

    // time vectors
    dVector _Cp;
    dVector _Cr;
    dVector _Cf;
    dVector _Cb;
    dVector _Ct;

    // DownSampled vectors
    dVector _CpBinned;
    dVector _CrBinned;
    dVector _CtBinned;

    void initializeBins();
    dVector downSample(dVector original);
    void addNoise(double noiseScale, dVector &downSampled);
    double GaussianRandomizer(double sigma, double cutoff);

    // calculate analytical solutions
    double integralOf(dVector tissue, int iTime);
    dVector differentiateTissueVector();
    dVector calculateConvolution(dVector tissue);

    void fitReferenceRegion();
    double baselineBasisFunction(int iPoly, double x);
    double gammaVariateFunction(double time, double onset, double alpha, double tau);


public:
    simEngine();
    void updateFineSamples();
    void run();
    void generatePlasmaTAC();
    void generateReferenceTAC();
    void generateTargetTAC();
    void generateTargetSRTM();
    void generateTargetFRTM();
    dVector FRTMOldAnalyticalSolution();
    dVector FRTMNewAnalyticalSolution();

    // setters
    inline void setSRTM()                     {_SRTM = true;}
    inline void setFRTM()                     {_SRTM = false;}
    inline void setNumberBins(int nBins)      {_nBins = nBins; initializeBins();}
    inline void setMagBolus(double value)     {_magBolus = value;}
    inline void setTauBolus(double value)     {_tauBolus = value;}
    inline void setKBol(double value)         {_KBol = value;}
    inline void setFastFraction(double value) {_fracFast = value;}
    inline void setKFastElim(double value)    {_kFast = value;}
    inline void setKSlowElim(double value)    {_kSlow = value;}
    inline void setTauFastElim(double value)  {_kFast = 1./value;}
    inline void setTauSlowElim(double value)  {_kSlow = 1./value;}
    inline void setK1Ref(double value)        {_K1Ref = value;}
    inline void setk2Ref(double value)        {_k2Ref = value;}
    inline void setTau1Ref(double value)      {_K1Ref = 1./value;}
    inline void setTau2Ref(double value)      {_k2Ref = 1./value;}
    inline void setk4(double value)           {_k4 = value;}
    inline void setTau4(double value)         {_k4 = 1./value;}
    inline void setR1(double value)           {_R1 = value;}
    inline void setBP0(double value)          {_BP0 = value;}
    inline void setK1(double value)           {_K1 = value;}
    inline void setk2(double value)           {_k2 = value;}
    inline void setk3(double value)           {_k3 = value;}
    inline void setNoiseRef(double value)     {_noiseRef = value;}
    inline void setNoiseTar(double value)     {_noiseTar = value;}
    inline void setChallengeTime(double value){_challengeTime  = value;}
    inline void setChallengeMag(double value) {_deltaBPPercent = value;}
    inline void setPlasmaPercentRef(double value){_percentPlasmaRef = value;}
    inline void setPlasmaPercentTar(double value){_percentPlasmaTar = value;}
    inline void setSamplesPerBin(int lBin, int nSamples)  {_numberSamplesPerBin[lBin] = nSamples; updateFineSamples();}
    inline void setDurationBin(int lBin, int duration) {_durationBinSec[lBin] = duration;         updateFineSamples();}
    inline void setDurationBins(iVector durationVector)   {_durationBinSec = durationVector;      updateFineSamples();}
    inline void setStartWithPlasma() {_startWithPlasma=true;}
    inline void setStartWithRR()     {_startWithPlasma=false;}

    // getters
    inline int getNumberTimeBinsFine()   {return _dtFine.size();}
    inline int getNumberTimeBinsCoarse() {return _timeCoarse.size();}
    inline bool getSRTM()           {return _SRTM;}
    inline int getNumberBins()      {return _nBins;}
    inline double getMagBolus()     {return _magBolus;}
    inline double getTauBolus()     {return _tauBolus;}
    inline double getKBol()  {return _KBol;}
    inline double getFastFraction() {return _fracFast;}
    inline double getKFastElim()    {return _kFast;}
    inline double getKSlowElim()    {return _kSlow;}
    inline double getTauFastElim()  {return 1./_kFast;}
    inline double getTauSlowElim()  {return 1./_kSlow;}
    inline double getK1Ref()        {return _K1Ref;}
    inline double getk2Ref()        {return _k2Ref;}
    inline double getTau1Ref()      {return 1./_K1Ref;}
    inline double getTau2Ref()      {return 1./_k2Ref;}
    inline double getk4()           {return _k4;}
    inline double getTau4()         {return 1./_k4;}
    inline double getR1()           {return _R1;}
    inline double getBP0()          {return _BP0;}
    inline double getK1()           {return _K1;}
    inline double getk2()           {return _k2;}
    inline double getk3()           {return _k3;}
    inline double getk2a()          {return _k2/(1.+_BP0);}
    inline double getk2k3()         {return _k2 * _k3;}  // or k2 * k4 * BPnd
    inline int getSamplesPerBin(int lBin) {return _numberSamplesPerBin[lBin];}
    inline double getDurationPerBin(int lBin) {return static_cast<double>(_durationBinSec[lBin])*60.;}
    inline double getDurationPerBinSec(int lBin) {return _durationBinSec[lBin];}
    inline iVector getTimeBinVectorSec()    {return _durationBinSec;}
    inline double getdtFine(int iTime) {return _dtFine[iTime];}
    inline double getTimeFine(int iTime) {return _timeFine[iTime];}
    inline double getTimeCoarse(int iTime) {return _timeCoarse[iTime];}
    inline double getdtCoarse(int iTime) {return _dtCoarse[iTime];}
    double getDurationScan();
    double getdk2a();
    double getdk2k3();

    inline double getNoiseRef()     {return _noiseRef;}
    inline double getNoiseTar()     {return _noiseTar;}
    inline double getChallengeTime(){return _challengeTime;}
    inline double getChallengeMag() {return _deltaBPPercent;}
    inline double getPlasmaPercentRef(){return _percentPlasmaRef;}
    inline double getPlasmaPercentTar(){return _percentPlasmaTar;}
    // concentations
    inline double getCp(int iTime)      {return _Cp[iTime];}
    inline double getCr(int iTime)      {return _Cr[iTime];}
    inline double getCt(int iTime)      {return _Ct[iTime];}
    inline double getCpDown(int iTime)  {return _CpBinned[iTime];}
    inline double getCrDown(int iTime)  {return _CrBinned[iTime];}
    inline double getCtDown(int iTime)  {return _CtBinned[iTime];}

};

#endif // SIMENGINE_H
