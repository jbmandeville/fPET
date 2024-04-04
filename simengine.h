#ifndef SIMENGINE_H
#define SIMENGINE_H

#include <QDebug>

#include "io.h"
#include "generalglm.h"

enum simulationStartingPoint
{
    simStart_fromPlasma,
    simStart_fromDataFit
};

class simEngine
{
private:
    // Input quantities ///////////////////////
    // Setup bins
    int _nBins = 90;
    iVector _durationBinSec;     // [_nBins]; seconds
    iVector _numberSamplesPerBin;  // [_nBins]
    dVector _dtFine;          // [sum of _numberSamplesPerBin = # of fine bins]; units = sec
    dVector _timeFine;        // [sum of _numberSamplesPerBin = # of fine bins]
    dVector _dtCoarse;        // [_nBins]
    dVector _timeCoarse;      // [_nBins]
    // starting point: plasma or RR
    // Infusion
    double _magBolus=150.;
    double _tauBolus=0.25;
    double _KBol=0.1;
    double _KBolDelay=0.;
    // Elimination
    double _fracFast=0.9;   // fraction of elimination that is fast
    double _kFast=1./1.2;    // 1/T_alpha fast elimination
    double _kSlow=1./20.;    // 1/T_beta  slow elimination
    double _percentPlasmaTar=0.0 ; // CBV fraction (set to 0 to ignore plasma)
    // Target region: state
    // Challenge
    double _challengeOnsetTime1=20.;
    double _challengeOffsetTime1=30.;
    double _challengeOnsetTime2=50.;
    double _challengeOffsetTime2=60.;

    // Where to start
    // structure for fitting the RR
    LOESS _quadLOESS;

    // derived quantities //////////////////////////
    double _K1=0.2;
    double _k3=0.2;    // k3 = BP * koff
    double _k2_div_k3=1.;
    double _dk3_percent1=30.;
    double _dk3_percent2=30.;
    // noise
    double _noiseTar=0.;

    // time vectors
    dVector _Cp;
    dVector _Cf;
    dVector _Cb;
    dVector _Ct;

    // DownSampled vectors
    dVector _CpBinned;
    dVector _CtBinned;

    void initializeBins();
    dVector downSample(dVector original);
    void addNoise(double noiseScale, dVector &downSampled);
    double GaussianRandomizer(double sigma, double cutoff);

public:
    simEngine();
    void updateFineSamples();
    void run();
    void generatePlasmaTAC();
    void generateTargetTAC();

    // setters
    inline void setNumberBins(int nBins)      {_nBins = nBins; initializeBins();}
    inline void setMagBolus(double value)     {_magBolus = value;}
    inline void setTauBolus(double value)     {_tauBolus = value;}
    inline void setKBol(double value)         {_KBol = value;}
    inline void setKBolDelay(double value)    {_KBolDelay = value;}
    inline void setFastFraction(double value) {_fracFast = value;}
    inline void setKFastElim(double value)    {_kFast = value;}
    inline void setKSlowElim(double value)    {_kSlow = value;}
    inline void setTauFastElim(double value)  {_kFast = 1./value;}
    inline void setTauSlowElim(double value)  {_kSlow = 1./value;}
    inline void setK1(double value)           {_K1 = value;}
    inline void setk2 (double value)          {_k2_div_k3 = value / _k3;}
    inline void setk3(double value)           {_k3 = value;}
    inline void setNoiseTar(double value)     {_noiseTar = value;}
    inline void setChallengeOnsetTime1(double value) {_challengeOnsetTime1  = value;}
    inline void setChallengeOffsetTime1(double value){_challengeOffsetTime1 = value;}
    inline void setChallengeOnsetTime2(double value) {_challengeOnsetTime2  = value;}
    inline void setChallengeOffsetTime2(double value){_challengeOffsetTime2 = value;}
    inline void setChallengeMag1(double value) {_dk3_percent1 = value;}
    inline void setChallengeMag2(double value) {_dk3_percent2 = value;}
    inline void setPlasmaPercentTar(double value){_percentPlasmaTar = value;}
    inline void setSamplesPerBin(int lBin, int nSamples)  {_numberSamplesPerBin[lBin] = nSamples; updateFineSamples();}
    inline void setDurationBin(int lBin, int duration) {_durationBinSec[lBin] = duration;         updateFineSamples();}
    inline void setDurationBins(iVector durationVector)   {_durationBinSec = durationVector;      updateFineSamples();}

    // getters
    inline int getNumberTimeBinsFine()   {return _dtFine.size();}
    inline int getNumberTimeBinsCoarse() {return _timeCoarse.size();}
    inline int getNumberBins()      {return _nBins;}
    inline double getMagBolus()     {return _magBolus;}
    inline double getTauBolus()     {return _tauBolus;}
    inline double getKBol()         {return _KBol;}
    inline double getKBolDelay()    {return _KBolDelay;}
    inline double getFastFraction() {return _fracFast;}
    inline double getKFastElim()    {return _kFast;}
    inline double getKSlowElim()    {return _kSlow;}
    inline double getTauFastElim()  {return 1./_kFast;}
    inline double getTauSlowElim()  {return 1./_kSlow;}
    inline double getK1()           {return _K1;}
    inline double getk2()           {return _k2_div_k3 * _k3;}
    inline double getk3()           {return _k3;}
    inline int getSamplesPerBin(int lBin) {return _numberSamplesPerBin[lBin];}
    inline double getDurationPerBin(int lBin) {return static_cast<double>(_durationBinSec[lBin])*60.;}
    inline double getDurationPerBinSec(int lBin) {return _durationBinSec[lBin];}
    inline iVector getTimeBinVectorSec()    {return _durationBinSec;}
    inline double getdtFine(int iTime) {return _dtFine[iTime];}
    inline double getTimeFine(int iTime) {return _timeFine[iTime];}
    inline double getTimeCoarse(int iTime) {return _timeCoarse[iTime];}
    inline double getdtCoarse(int iTime) {return _dtCoarse[iTime];}

    inline double getNoiseTar()     {return _noiseTar;}
    inline double getChallengeOnsetTime1()  {return _challengeOnsetTime1;}
    inline double getChallengeOffsetTime1() {return _challengeOffsetTime1;}
    inline double getChallengeOnsetTime2()  {return _challengeOnsetTime2;}
    inline double getChallengeOffsetTime2() {return _challengeOffsetTime2;}
    inline double getChallengeMagPercent1() {return _dk3_percent1;}
    inline double getChallengeMagPercent2() {return _dk3_percent2;}
    inline double getPlasmaPercentTar(){return _percentPlasmaTar;}
    // concentations
    inline double getCp(int iTime)      {return _Cp[iTime];}
    inline double getCt(int iTime)      {return _Ct[iTime];}
    inline double getCpCoarse(int iTime)  {return _CpBinned[iTime];}
    inline double getCtCoarse(int iTime)  {return _CtBinned[iTime];}
    inline dVector getCtCoarse() {return _CtBinned;}
    inline dVector getTimeCourse()   {return _timeCoarse;}
};

#endif // SIMENGINE_H
