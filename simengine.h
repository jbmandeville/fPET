#ifndef SIMENGINE_H
#define SIMENGINE_H

#include "io.h"


class simEngine
{
private:
    // Input quantities ///////////////////////
    // Setup
    double _duration=90.;
    double _dt=0.1;
    int _lDownSample=10;
    // Infusion
    double _magBolus=1.e5;
    double _tauBolus=1.5;
    double _alphaBolus=1.0;
    double _magInfusion=1.e3;
    // Elimination
    double _fracFast=0.95;   // fraction of elimination that is alpha
    double _kFast=1./1.5;    // 1/T_alpha fast elimination
    double _kSlow=1.;        // 1/T_beta  slow elimination
    double _fracPlasma=0.0 ; // CBV fraction (set to 0 to ignore plasma)
    // Reference region
    double _K1Ref=0.15;
    double _k2Ref=1./2.;
    // Target region: state
    double _k4=1./10.;
    double _R1=1.;
    double _BP0=4.;     // BP = k3/koff = kon*Bavail0/koff
    // Challenge
    double _challengeTime=40.;
    double _deltaBPPercent=10.;  // deltaBPnd in %

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

    dVector downSample(dVector original);
    void addNoise(double noiseScale, dVector &downSampled);
    double GaussianRandomizer(double sigma, double cutoff);

public:
    simEngine();
    void run();
    void generatePlasmaTAC();
    void generateReferenceTAC();
    void generateTargetTAC();

    // setters
    inline void setDuration(double value)     {_duration = value;}
    inline void setStepSize(double value)     {_dt = value;}
    inline void setDownSampling(int lValue)   {_lDownSample = lValue;}
    inline void setMagBolus(double value)     {_magBolus = value;}
    inline void setTauBolus(double value)     {_tauBolus = value;}
    inline void setAlphaBolus(double value)   {_alphaBolus = value;}
    inline void setMagInfusion(double value)  {_magInfusion = value;}
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

    // getters
    inline double getDuration()     {return _duration;}
    inline double getStepSize()     {return _dt;}
    inline int getDownSampling()    {return _lDownSample;}
    inline double getMagBolus()     {return _magBolus;}
    inline double getTauBolus()     {return _tauBolus;}
    inline double getAlphaBolus()   {return _alphaBolus;}
    inline double getMagInfusion()  {return _magInfusion;}
    inline double getFastFraction() {return _fracFast;}
    inline double getKFastElim()    {return _kFast;}
    inline double getKSlowElim()    {return _kSlow;}
    inline double getTauFastElim()  {return 1./_kFast;}
    inline double getTauSlowElim()  {return 1./_kFast;}
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
    inline double getNoiseRef()     {return _noiseRef;}
    inline double getNoiseTar()     {return _noiseTar;}
    inline double getChallengeTime(){return _challengeTime;}
    inline double getChallengeMag() {return _deltaBPPercent;}
    // concentations
    inline double getCp(int iTime)      {return _Cp[iTime];}
    inline double getCr(int iTime)      {return _Cr[iTime];}
    inline double getCt(int iTime)      {return _Ct[iTime];}
    inline double getCpDown(int iTime)  {return _CpBinned[iTime];}
    inline double getCrDown(int iTime)  {return _CrBinned[iTime];}
    inline double getCtDown(int iTime)  {return _CtBinned[iTime];}

};

#endif // SIMENGINE_H
